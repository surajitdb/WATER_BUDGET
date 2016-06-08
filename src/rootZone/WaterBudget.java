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

package rootZone;

import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static org.jgrasstools.gears.libs.modules.JGTConstants.isNovalue;

import java.util.HashMap;
import java.util.Set;
import java.util.Map.Entry;

import oms3.annotations.Description;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Out;
import oms3.annotations.Unit;

import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.feature.SchemaException;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.jgrasstools.gears.utils.coverage.CoverageUtilities;
import org.jgrasstools.hortonmachine.modules.statistics.cb.OmsCb;

import java.awt.image.RenderedImage;
import java.awt.image.WritableRaster;
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


	@Description("Input rain Hashmap")
	@In
	public HashMap<Integer, double[]> inRainValues;


	@Description("Input rain Hashmap")
	@In
	public HashMap<Integer, double[]> inSnowValues;

	@Description("ET: ET value for the given time considered")
	double ET;

	@Description("Input ET Hashmap")
	@In
	public HashMap<Integer, double[]> inETvalues;

	@Description("Input Discharge Hashmap: first contribution from other HRUs")
	@In
	public HashMap<Integer, double[]> inDischargevalues1;

	@Description("Input Discharge Hashmap: second contribution from other HRUs")
	@In
	public HashMap<Integer, double[]> inDischargevalues2;


	@Description("time step of the simulation")
	@In
	public int inTimestep;


	@Description("Integration time")
	double dt ;

	@Description("ID of the basin")
	@In
	public static int ID ;



	@Description("Parameter of the non-linear Reservoir model "
			+ "for the considered layer")
	@In
	public static double a ;


	@Description("Parameter of non-linear reservoir, for the upper layer")
	@In
	public static double c;


	@Description("Pore volume in the root zone")
	@In
	public static double nZ;

	@Description("partitioning coefficient between the reserovir")
	@Unit("-")
	@In
	@Out
	public double alpha;


	@Description("Degree of spatial variability of the soil moisture capacity.")
	@In
	@Unit("-")
	public Double pB;

	@Description("Maximum value of the water storage, needed for the"
			+ "computation of the Actual EvapoTraspiration")
	@In
	public static double s_max;


	@Description("Maximum recharge rate of the lower layer")
	@In
	public static double Re;


	@Description("Discharge model: NonLinearReservoir, Clapp-H")
	@In
	public String Q_model;

	@Description("ET model: AET,ExternalValues, nullET")
	@In
	public String ET_model;

	@Description("ODE solver ")
	@In
	public String solver_model;

	@Description("Simluted value of AET"
			+ "at a given time step")
	double AET;

	DischargeModel Qmodel;
	ETModel ETmodel;

	@Description("Rescaled distance map")
	@In
	public GridCoverage2D inRescaledsup = null;

	@Description("Channel celerity")
	@Unit("m/s")
	@In
	public double pCelerity;

	@Description("The output HashMap with the Water Storage  ")
	@Out
	public HashMap<Integer, double[]> outHMStorage= new HashMap<Integer, double[]>() ;

	@Description("The output HashMap with the discharge ")
	@Out
	public HashMap<Integer, double[]> outHMDischarge= new HashMap<Integer, double[]>() ;

	@Description("The output HashMap with the AET ")
	@Out
	public HashMap<Integer, double[]> outHMEvapotranspiration = new HashMap<Integer, double[]>() ;


	@Description("The output HashMap with the outflow "
			+ "which drains to the lower layer")
	@Out
	public HashMap<Integer, double[]> outHMR= new HashMap<Integer, double[]>() ;

	double S_i;
	int step;



	/**
	 * Process: reading of the data, computation of the
	 * storage and outflows
	 *
	 * @throws Exception the exception
	 */
	@Execute
	public void process() throws Exception {
		checkNull(inRainValues);

		/**Input data reading*/
		double rain = inRainValues.get(ID)[0];
		if (isNovalue(rain)) rain= 0;

		double snow = 0;					
		if (inSnowValues != null) snow=inSnowValues.get(ID)[0];
		if (isNovalue(snow)) snow= 0;

		double Qinput1=0;
		if (inDischargevalues1 != null) Qinput1 =inDischargevalues1.get(ID)[0];
		if (isNovalue(Qinput1)) Qinput1= 0;

		double Qinput2=0;
		if (inDischargevalues2 != null) Qinput2 =inDischargevalues2.get(ID)[0];
		if (isNovalue(Qinput2)) Qinput2= 0;


		double totalInputFluxes=rain+snow+Qinput1+Qinput2;

		ET=0;
		if (inETvalues != null) ET = inETvalues.get(ID)[0];
		if (isNovalue(ET)) ET= 0;

		double alpha=alpha(S_i,totalInputFluxes);
		double waterStorage=computeS((1-alpha)*totalInputFluxes);
		double discharge=computeQ(waterStorage);
		double evapotranspiration=computeAET(S_i);
		double drainage=computeR(discharge);


		/** Save the result in  hashmaps for each station*/
		storeResult_series(ID,waterStorage,discharge,evapotranspiration,drainage);


		step++;

	}

	/**
	 * Compute the water storage
	 *
	 * @param J: input rain 
	 * @param Qinput : input discharge
	 * @param ET: input potential ET
	 * @return the water storage, according to the model and the layer
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	public double computeS(double totalInputFluxes) throws IOException {
		/**integration time*/
		dt=1E-4;

		/** SimpleFactory for the computation of Q, according to the model*/
		Qmodel=SimpleDischargeModelFactory.createModel(Q_model, a, S_i, c);
		double Qmod=Qmodel.dischargeValues();

		/** SimpleFactory for the computation of ET, according to the model*/
		ETmodel=SimpleETModelFactory.createModel(ET_model,ET,S_i,s_max);
		double ETmod=ETmodel.ETValues();

		/** Creation of the differential equation*/
		FirstOrderDifferentialEquations ode=new waterBudgetODE(nZ,totalInputFluxes,Qmod, ETmod);			

		/** Boundaries conditions*/
		double[] y = new double[] { S_i, 0 };

		/** Choice of the ODE solver */	
		SolverODE solver;
		solver=SimpleIntegratorFactory.createSolver(solver_model, dt, ode, y);

		/** result of the resolution of the ODE, if nZ=1, S_i=S_t
		 * and setting of the new initial condition (S_i)*/
		S_i=solver.integrateValues();
		double S_t=S_i*nZ;

		/** Check of the Storage values: they cannot be negative*/
		if (S_t<0) S_t=0;

		return S_t;
	}

	/**
	 * Compute computation of the discharge according to the mode:
	 * mode external --> external value
	 * else --> non-linear reservoir model
	 *
	 * @param Qinput: input discharge value
	 * @return the double value of the simulated discharge
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	public double computeQ( double S_i) throws IOException {
		Qmodel=SimpleDischargeModelFactory.createModel(Q_model, a, S_i, c);
		double Q=Qmodel.dischargeValues();
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
		double R=Math.min(Q, Re);
		return R;
	}

	/**
	 * Compute the AET
	 *
	 * @param ETinput: the input potential ET
	 * @return the double value od the AET
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	public double computeAET(double S_i) throws IOException {
		ETmodel=SimpleETModelFactory.createModel(ET_model,ET,S_i,s_max);
		AET=ETmodel.ETValues();
		return AET;
	}


	private double alpha( double S_i_prev, double Pval) {
		double pCmax=s_max *(pB+1);
		double coeff1 = ((1.0 - ((pB + 1.0) * (S_i_prev) / pCmax)));
		double exp = 1.0 / (pB + 1.0);
		double ct_prev = pCmax * (1.0 - pow(coeff1, exp));
		double UT1 = max((Pval - pCmax + ct_prev), 0.0);     
		return alpha=UT1/Pval;
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
			double evapotranspiration,double drainage) throws SchemaException {

		outHMStorage.put(ID, new double[]{waterStorage});
		outHMDischarge.put(ID, new double[]{discharge});
		outHMEvapotranspiration.put(ID, new double[]{evapotranspiration});
		outHMR.put(ID, new double[]{drainage});

	}


}