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
 * @description The Class WaterBudget solves the water budget equation,
 *              according to the models chosen for the simulation of the
 *              discharge and of the AET. The outputs of the class are 4
 *              hashmaps with the water storage values, simulated discharge
 *              values, simulated AET and the quick discharge (if there is any).
 *              The simulated discharge is partitioned in two flows: a flow
 *              which drains to the lower layer and the quick discharge. For the
 *              lower layer, we don't have nor ET processes and the partitioning
 *              of the discharge in quick and drainage to other layers, so the
 *              third and fourth columns of the first hashmap and the second
 *              hashmap will be equal to zero.
 *
 * @todo Check the name of the variables, defining the proper name where
 *       required
 * @author Marialaura Bancheri
 * @author sidereus (francesco.serafin.3@gmail.com)
 * @version
 */
public class WaterBudget extends JGTModel{


	@Description("precipitation: precipitation value for the given time considered")
	double precipitation;

	@Description("Input Precipitation Hashmap")
	@In
	public HashMap<Integer, double[]> inPrecipValues;

	@Description("ET: ET value for the given time considered")
	double evapotranspiration;

	@Description("Input ET Hashmap")
	@In
	public HashMap<Integer, double[]> inEvapotranspValues;

	@Description("Input Discharge Hashmap: this is the measured "
			+ "discharge time series, in the case of resultion of the "
			+ "water budget with given the external values")
	@In
	public HashMap<Integer, double[]> inDischargeValues;

	@Description("time step of the simulation")
	@In
	public int inTimestep;

	@Description("Integration time")
	double dt; // proposed name: integrationTime

	@Description("Area of the basin")
	@In
	public double basinArea;

	@Description("Parameter of the non-linear Reservoir model "
			+ "for the considered layer")
	@In
	public double a; // define a proper name for variable a

	@Description("Parameter of non-linear reservoir, for the upper layer")
	@In
	public double b; // define a proper name for variable b

	@Description("Pore volume in the root zone")
	@In
	public double poreVolumeInRootZone;

	@Description("Maximum value of the water storage, needed for the"
			+ "computation of the Actual EvapoTraspiration")
	@In
	public double waterStorageMaxValue;

	@Description("Maximum recharge rate of the lower layer")
	@In
	public double Re; // proposed name: rechargeRateMaxValueLowerLayer

	@Description("Discharge model: NonLinearReservoir, ExternalValues")
	@In
	public String dischargeModelName;

	@Description("ET model: AET,ExternalValues")
	@In
	public String evapotranspirationModelName;

	@Description("ODE solver ")
	@In
	public String odeSolverModelName;

	@Description("water storage initial condition")
	@In
	public static double waterStorageInitalConditions;

	@Description("Simluted value of Water storage,"
			+ "at a given time step")
	double S_t; // proposed name: waterStorage

	@Description("Outflow from the upper layer, which drains to"
			+ "the lower layer")
	double R; // proposed name: flowUpperToLowerLayer

	@Description("Simluted value of quick flow"
			+ "at a given time step")
	double quickFlow;

	@Description("Simluted value of AET"
			+ "at a given time step")
	double actualEvapotranspiration;

	DischargeModel dischargeModel;
	EvapotranspirationModel evapotranspirationModel;

	@Description("The output HashMap with the Water Storage")
	@Out
	public HashMap<Integer, double[]> outHMStorage; // proposed name: outWaterStorage

	@Description("The output HashMap with the discharge ")
	@Out
	public HashMap<Integer, double[]> outHMDischarge; // proposed name: outDischarge

	@Description("The output HashMap with the AET ")
	@Out
	public HashMap<Integer, double[]> outHMEvapotranspiration; // proposed name: outEvapotranspiration

	@Description("The output HashMap with the  quick runoff")
	@Out
	public HashMap<Integer, double[]> outHMQuick; // proposed name: outQuickDischarge

	@Description("The output HashMap with the outflow "
			+ "which drains to the lower layer")
	@Out
	public HashMap<Integer, double[]> outHMR; // proposed name: outFlowUpperToLowerLayer

	/**
	 * @brief Default constructor
	 */
	public WaterBudget() {}

	/**
	 * @description Process: reading of the data, computation of the storage and
	 *              outflows
	 *
	 * @throws Exception
	 *             the exception
	 * @todo Better implementation of the -Input data reading- section. Not
	 *       clear why just the first value is retrieved
	 */
	@Execute
	public void process() throws Exception {
		checkNull(inPrecipValues);

		// reading the ID of all the stations 
		Set<Entry<Integer, double[]>> entrySet = inPrecipValues.entrySet();

		allocationOutputVariables(entrySet.size());

		// iterate over the station
		for( Entry<Integer, double[]> entry : entrySet ) {
			Integer ID = entry.getKey();

			/**Input data reading*/
			precipitation = inPrecipValues.get(ID)[0];
			if (isNovalue(precipitation)) precipitation= 0;

			double Qinput=0;
			if (inDischargeValues != null) Qinput =inDischargeValues.get(ID)[0];


			evapotranspiration=0;
			if (inEvapotranspValues != null) evapotranspiration = inEvapotranspValues.get(ID)[0];



			computeWaterStorage(Qinput);
			double discharge = computeDischarge(Qinput);
			computeActualEvapotranspiration();
			computeQuickRunoff(discharge);
			computeR(discharge); // proposed name for method: computeFlowTowardLayers

			/** Save the result in  hashmaps for each station*/
			storeResult_series(ID,S_t,discharge,actualEvapotranspiration,quickFlow,R);

		}

	}

	/**
	 * Compute the water storage
	 *
	 * @param Qinput : input discharge
	 * @return the water storage, according to the model and the layer
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	public void computeWaterStorage(double Qinput) throws IOException {
		/**integration time*/
		dt=1E-4;

		/** SimpleFactory for the computation of Q, according to the model*/
		dischargeModel=SimpleDischargeModelFactory.createModel(dischargeModelName, Qinput, basinArea, a, waterStorageInitalConditions, b);
		double Qmod=dischargeModel.dischargeValues();

		/** SimpleFactory for the computation of ET, according to the model*/
		evapotranspirationModel=SimpleETModelFactory.createModel(evapotranspirationModelName,evapotranspiration,waterStorageInitalConditions,waterStorageMaxValue);
		double ETmod=evapotranspirationModel.ETValues();

		/** Creation of the differential equation*/
		FirstOrderDifferentialEquations ode=new waterBudgetODE(poreVolumeInRootZone,precipitation, Qmod, ETmod);

		/** Boundaries conditions*/
		double[] y = new double[] { waterStorageInitalConditions, 0 };

		/** Choice of the ODE solver */	
		SolverODE solver;
		solver=SimpleIntegratorFactory.createSolver(odeSolverModelName, dt, ode, y);

		/** result of the resolution of the ODE, if nZ=1, S_i=S_t
		 * and setting of the new initial condition (S_i)*/
		waterStorageInitalConditions=solver.integrateValues();
		S_t=waterStorageInitalConditions*poreVolumeInRootZone;

		/** Check of the Storage values: they cannot be negative*/
		if (S_t<0) S_t=0;
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
	public double computeDischarge( double Qinput) throws IOException {
		dischargeModel=SimpleDischargeModelFactory.createModel(dischargeModelName, Qinput, basinArea, a, waterStorageInitalConditions, b);
		return dischargeModel.dischargeValues();
	}

	/**
	 * Compute the outflow toward the lower layer
	 *
	 * @param Q: simulated discharge for the considered layer
	 * @return the double value of the outflow toward the lower layer
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	public void computeR(double Q) throws IOException {
		R = Math.min(Q, Re);
	}

	/**
	 * Compute quick runoff from the layer.
	 *
	 * @param Q: simulated discharge for the considered layer
	 * @return the double value of the quick runoff from the layer
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	public void computeQuickRunoff(double Q) throws IOException {
		quickFlow = Q-R;
	}

	/**
	 * Compute the AET
	 *
	 * @param ETinput: the input potential ET
	 * @return the double value od the AET
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	public void computeActualEvapotranspiration() throws IOException {
		evapotranspirationModel=SimpleETModelFactory.createModel(evapotranspirationModelName,evapotranspiration,waterStorageInitalConditions,waterStorageMaxValue);
		actualEvapotranspiration = evapotranspirationModel.ETValues();
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

	private void allocationOutputVariables(final int size) {

		outHMStorage = new HashMap<Integer, double[]>(size);
		outHMDischarge = new HashMap<Integer, double[]>(size);
		outHMEvapotranspiration = new HashMap<Integer, double[]>(size);
		outHMQuick = new HashMap<Integer, double[]>(size);
		outHMR = new HashMap<Integer, double[]>(size);

	}

}
