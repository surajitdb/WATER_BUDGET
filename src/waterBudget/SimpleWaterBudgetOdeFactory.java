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

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations;

// TODO: Auto-generated Javadoc
/**
 * A simple factory pattern that create concrete ODE
 * for the resolution of the water budget equations: FirstLayer, 
 * SecondLayer and ExternalValues. 
 * @author Marialaura Bancheri
 */
public class SimpleWaterBudgetOdeFactory {

	
	/**
	 * Creates three WaterBudget Ode.
	 *
	 * @param type: the chosen layer 
	 * @param a: the coefficient of the non-linear reservoir
	 * @param b: the exponent of the non-linear reservoir
	 * @param nZ: the product of the porosity and the depth of the root zone
	 * @param s_max: the maximum value of the soil moisture
	 * @param A : the area of the basin
	 * @param J: the precipitation value
	 * @param ET : the ET values
	 * @param Q: the discharge value
	 * @param S: the soil moisture/water storage at the previous time-step
	 * @return the first order differential equations
	 */
	public static FirstOrderDifferentialEquations createOde(String type,double Qmod,double nZ,double ETmod, 
			double J){
		FirstOrderDifferentialEquations ode=null;
		
		if (type.equals("WB")){
			ode=new Layer(nZ, J, Qmod,ETmod);
		}
			
		return ode;
		
	}

}
