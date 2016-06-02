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

package sumSeries;

import static org.jgrasstools.gears.libs.modules.JGTConstants.isNovalue;

import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Set;
import java.util.Map.Entry;

import oms3.annotations.Description;
import oms3.annotations.Execute;
import oms3.annotations.In;
import oms3.annotations.Out;

import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.feature.FeatureIterator;
import org.geotools.feature.SchemaException;
import org.jgrasstools.gears.libs.modules.JGTModel;
import org.opengis.feature.simple.SimpleFeature;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;
import org.jgrasstools.hortonmachine.modules.hydrogeomorphology.adige.core.PfafstetterNumber;

import java.io.IOException;




// TODO: Auto-generated Javadoc
/**
 * The Class SumSeries compute the sum of two time series
 */
public class SumSeriesConnections extends JGTModel{

	@Description("Connetction between HRUS")
	@In
	public int [] connections;

	
	@Description("Input upper layer discharge Hashmap")
	@In
	public HashMap<Integer, double[]> inQ;
	

	@Description("Sum at a given time step")
	double Qtot;
	
	int ID ;
	double sum;

	@Description("The output HashMap with the sum"
			+ "for the considered layer ")
	@Out
	public HashMap<Integer, double[]> outHMQtot= new HashMap<Integer, double[]>() ;



	/**
	 * Process.
	 *
	 * @throws Exception the exception
	 */
	@Execute
	public void process() throws Exception {
		
		checkNull(inQ);
		


		for (int i=0; i<connections.length;i++){
							
			ID=connections[i];

			double Q2 =inQ.get(ID)[0];
			if (isNovalue(Q2)) Q2= 0;

			/** sum of the given quantities */
			sum=sum(sum,Q2);

		}
		
		/** Save the result in hashmap output */
		storeResult_series(ID,sum);
		sum=0;

		}




	/**
	 * Sum:	computation of the sum
	 *
	 * @param Q1: the first quantity
	 * @param Q2 : the second quantity
	 * @return the double value of the sum
	 * @throws IOException Signals that an I/O exception has occurred.
	 */
	public double sum(double Q1, double Q2) throws IOException {
		Qtot=Q1+Q2;
		return Qtot;
	}


	/**
	 * Store result_series
	 *
	 * @param Qtot: the vector with the sums
	 * @throws SchemaException the schema exception
	 */
	private void storeResult_series(int ID,double Qtot) throws SchemaException {	
		outHMQtot.put(ID, new double[]{Qtot});


	}
}