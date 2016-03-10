/*
 * GNU GPL v3 License
 *
 * Copyright 2015 Marialaura Bancheri
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package waterBudget;

import static org.jgrasstools.gears.libs.modules.JGTConstants.isNovalue;

import java.util.HashMap;
import java.util.Set;
import java.util.Map.Entry;

import oms3.annotations.Description;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Out;

import org.geotools.feature.SchemaException;
import org.jgrasstools.gears.libs.modules.JGTModel;

import java.io.IOException;

import org.apache.commons.math3.ode.*;


/**
 * The Class WaterBudget solves the water budget equation, according to
 * the models chosen for the simulation of the discharge and of the AET. 
 * The outputs of the class are 4 hashmaps with the water storage values,
 * simulated discharge values, simulated AET and the quick discharge (if there
 * is any). The simulated discharge is partitioned in two flows: a flow which drains
 * to the lower layer and the quick discharge.
 * For the lower layer, we don't have nor ET processes and the 
 * partitioning of the discharge in quick and drainage to other layers, so the third and
 * fourth columns of the first hashmap and the second hashmap will be equal to zero.
 * @author Marialaura Bancheri
 */
public class WaterBudget extends JGTModel{


	@Description("precipitation: precipitation value for the given time considered") 
	double precipitation;

	@Description("Input Precipitation Hashmap")
	@In
	public HashMap<Integer, double[]> inPrecipvalues;

	@Description("ET: ET value for the given time considered")
	double ET;

	@Description("Input ET Hashmap")
	@In
	public HashMap<Integer, double[]> inETvalues;

	@Description("Input Discharge Hashmap: this is the measured "
			+ "discharge time series, in the case of resultion of the "
			+ "water budget with given the external values")
	@In
	public HashMap<Integer, double[]> inDischargevalues;


	@Description("time step of the simulation")
	@In
	public int inTimestep;


	@Description("Integration time")
	double dt ;


	@Description("Area of the basin")
	@In
	public static double A ;


	@Description("Parameter of the non-linear Reservoir model "
			+ "for the considered layer")
	@In
	public static double a ;


	@Description("Parameter of non-linear reservoir, for the upper layer")
	@In
	public static double b;


	@Description("Pore volume in the root zone")
	@In
	public static double nZ;


	@Description("Maximum value of the water storage, needed for the"
			+ "computation of the Actual EvapoTraspiration")
	@In
	public static double s_max;


	@Description("Maximum recharge rate of the lower layer")
	@In
	public static double Re;


	@Description("ODE mode: NonLinearReservoir, ExternalValues, nullQ")
	@In
	public String Q_model;

	@Description("ODE mode: AET,ExternalValues, nullET")
	@In
	public String ET_model;


	@Description("ODE solver ")
	@In
	public String solver_model;


	@Description("water storage initial condition ")
	@In
	public static double S_i;


	@Description("Simluted value of Water storage,"
			+ "at a given time step")
	double S_t;

	@Description("Outflow from the upper layer, which drains to"
			+ "the lower layer")
	double R; 

	@Description("Simluted value of quick flow"
			+ "at a given time step")
	double Qquick;

	@Description("Simluted value of AET"
			+ "at a given time step")
	double AET;

	DischargeModel model;
	ETModel ETmodel;

	@Description("The output HashMap with the Water Storage  ")
	@Out
	public HashMap<Integer, double[]> outHMStorage= new HashMap<Integer, double[]>() ;

	@Description("The output HashMap with the discharge ")
	@Out
	public HashMap<Integer, double[]> outHMDischarge= new HashMap<Integer, double[]>() ;

	@Description("The output HashMap with the AET ")
	@Out
	public HashMap<Integer, double[]> outHMEvapotranspiration = new HashMap<Integer, double[]>() ;

	@Description("The output HashMap with the  quick runoff")
	@Out
	public HashMap<Integer, double[]> outHMQuick= new HashMap<Integer, double[]>() ;

	@Description("The output HashMap with the outflow "
			+ "which drains to the lower layer")
	@Out
	public HashMap<Integer, double[]> outHMR= new HashMap<Integer, double[]>() ;


	/**
	 * Process: reading of the data, computation of the
	 * storage and outflows
	 *
	 * @throws Exception the exception
	 */
	@Execute
	public void process() throws Exception {
		checkNull(inPrecipvalues,inDischargevalues,inETvalues);

		// reading the ID of all the stations 
		Set<Entry<Integer, double[]>> entrySet = inPrecipvalues.entrySet();

		// iterate over the station
		for( Entry<Integer, double[]> entry : entrySet ) {
			Integer ID = entry.getKey();

			/**Input data reading*/
			precipitation = inPrecipvalues.get(ID)[0];
			if (isNovalue(precipitation)) precipitation= 0;

			double Qinput =inDischargevalues.get(ID)[0];
			if (isNovalue(Qinput)) Qinput= 0;

			if (s_max==0) ET=0;
			else { ET = inETvalues.get(ID)[0];
			if (isNovalue(ET)) ET= 0;}


			double waterStorage=computeS(Qinput);
			double discharge=computeQ(Qinput);
			double evapotranspiration=computeAET();	
			double quickRunoff=computeQuick(discharge);
			double drainage=computeR(discharge);

			/** Save the result in  hashmaps for each station*/
			storeResult_series(ID,waterStorage,discharge,evapotranspiration,quickRunoff,drainage);

		}




	}

	/**
	 * Compute the water storage
	 *
	 * @param J: input precipitation 
	 * @param Qinput : input discharge
	 * @param ET: input potential ET
	 * @return the water storage, according to the model and the layer
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	public double computeS(double Qinput) throws IOException {
		/**integration time*/
		dt=1E-4;

		/** SimpleFactory for the computation of Q, according to the model*/
		model=SimpleDischargeModelFactory.createModel(Q_model, Qinput, A, a, S_i, b);
		double Qmod=model.dischargeValues();

		/** SimpleFactory for the computation of ET, according to the model*/
		ETmodel=SimpleETModelFactory.createModel(ET_model,ET,S_i,s_max);
		double ETmod=ETmodel.ETValues();

		/** SimpleFactory for the creation of the differential equation*/
		FirstOrderDifferentialEquations ode;			
		ode=SimpleWaterBudgetOdeFactory.createOde("WB",Qmod,nZ,ETmod,precipitation);

		/** Boundaries conditions*/
		double[] y = new double[] { S_i, 0 };

		/** Choice of the ODE solver */	
		SolverODE solver;
		solver=SimpleIntegratorFactory.createSolver(solver_model, dt, ode, y);

		/** result of the resolution of the ODE, if nZ=1, S_i=S_t
		 * and setting of the new initial condition (S_i)*/
		S_i=solver.integrateValues();
		S_t=S_i*nZ;

		/** Check of the Storage values: they cannot be negative*/
		if (S_t<0) S_t=0;

		return S_t;
	}

	// computation of the discharge according to the mode: 
	// mode external --> external value
	// else --> model
	/**
	 * Compute computation of the discharge according to the mode:
	 * mode external --> external value
	 * else --> non-linear reservoir model
	 *
	 * @param Qinput: input discharge value
	 * @return the double value of the simulated discharge
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	public double computeQ( double Qinput) throws IOException {
		model=SimpleDischargeModelFactory.createModel(Q_model, Qinput, A, a, S_i, b);
		double Q=model.dischargeValues();
		return Q;
	}


	/**
	 * Compute the outflow toward the lower layer
	 *
	 * @param Q: simulated discharge for the considered layer
	 * @return the double value of the outflow toward the lower layer
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	public double computeR(double Q) throws IOException {
		R=Math.min(Q, Re);
		return R;
	}

	/**
	 * Compute quick runoff from the layer.
	 *
	 * @param Q: simulated discharge for the considered layer
	 * @return the double value of the quick runoff from the layer
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	public double computeQuick(double Q) throws IOException {
		Qquick=Q-R;
		return Qquick;
	}

	/**
	 * Compute the AET
	 *
	 * @param ETinput: the input potential ET
	 * @return the double value od the AET
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	public double computeAET() throws IOException {
		ETmodel=SimpleETModelFactory.createModel(ET_model,ET,S_i,s_max);
		AET=ETmodel.ETValues();
		return AET;
	}	


	/**
	 * Store of the results in hashmaps 
	 *
	 * @param waterStorage is the water storage
	 * @param discharge is the discharge
	 * @param evapotranspiration is the evapotranspiration
	 * @param quickRunoff is the water quick runoff from the layer
	 * @param drainage is drainage toward the lower layer
	 * @throws SchemaException the schema exception
	 */
	private void storeResult_series(int ID, double waterStorage,double discharge,
			double evapotranspiration, double quickRunoff,double drainage) throws SchemaException {

		outHMStorage.put(ID, new double[]{waterStorage});
		outHMDischarge.put(ID, new double[]{discharge});
		outHMEvapotranspiration.put(ID, new double[]{evapotranspiration});
		outHMQuick.put(ID, new double[]{quickRunoff});
		outHMR.put(ID, new double[]{drainage});


	}
}