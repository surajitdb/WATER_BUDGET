package runoff;


import org.jgrasstools.gears.libs.modules.ModelsEngine;


/**
 * The Class WFIUHKinematic.
 */
public class WFIUHKinematic {

	/** The width function */
	double [][] widthFunction;
	
	/** The total fluxes in input (rain+snow) */
	double totalInputFluxes;
	
	/** The input time step in minutes */
	int inTimestep;
	
	/** The area. */
	double area;
	
	/** The celerity of the channel */
	double pCelerity;
	
	/** The vector with discharge values at prevoius time step */
	double [] Q_i;
	
	/** The lenght of the input time series */
	int dim;
	
	/** The actual step. */
	int step;

	/**
	 * Instantiates a new IUH kinematic.
	 *
	 * @param widthFunction is the width function matrix
	 * @param totalInputFluxes are the totalInputFluxes for the time step
	 * @param inTimestep is the input time step in minutes
	 * @param area is the area
	 * @param pCelerity is the celerity of the channel
	 * @param Q_i the discharge at the prevoius time step
	 * @param dim is the length of the input time series
	 * @param step the actual computation step
	 */
	public WFIUHKinematic(double [][] widthFunction, double totalInputFluxes, int inTimestep, 
							double area, double pCelerity, double [] Q_i, int dim, int step){
		this.totalInputFluxes=totalInputFluxes;
		this.widthFunction=widthFunction;
		this.inTimestep=inTimestep;
		this.area=area;
		this.pCelerity=pCelerity;
		this.Q_i=Q_i;
		this.dim=dim;
		this.step=step;
	}

	/**
	 * Calculate the discharge time series.
	 *
	 * @return the double[]
	 */
	public double [] calculateQ() {

		// tcorr is the concentration time in seconds 
		double tcorr = widthFunction[widthFunction.length - 1][0];

		// is the duration of the precipitation in seconds
		int tpmax =inTimestep*60;

		//is the discharge computed in m^3/s according to the WFIUH at time step i+1 (i1)
		double[] Q_i1 = new double[(int)tcorr];


		for( int t = 1; t < tcorr-1; t += 1 ) {

			if (t <= tpmax) {

				Q_i1[t]=(double) (totalInputFluxes * area* ModelsEngine.width_interpolate(widthFunction, t, 0, 2))*pCelerity;

			} else {
				Q_i1[t]= (double) (totalInputFluxes * area* (ModelsEngine.width_interpolate(widthFunction, t, 0, 2) - ModelsEngine
						.width_interpolate(widthFunction, t - tpmax, 0, 2)))*pCelerity;
			}
		}

		// is the average discharge in m^3/s over 60 seconds 
		double [] Q=computeMean(Q_i1);

		// sum of the different contributes of the discharge (previous time step and actual time step)
		// where the two time series overlap
		for( int t = 0; t <  Q.length; t += 1 ) {	
			Q_i[t+step*inTimestep]= Q_i[t+step*inTimestep]+Q[t];
		}

		return Q_i;


	}


	/**
	 * Compute the mean of the discharge computed each second,
	 * in order to pass from seconds to minutes and have a faster code.
	 *
	 * @param runoff is the runoff computed each second with the WFIUH
	 * @return the double[] vector of the average runoff in a minute
	 */
	public double [] computeMean (double [] runoff){

		int step=60;

		
		double [] sum=new double [runoff.length/step];

		for(int t=0; t<sum.length;t++){
			for (int i=t*step;i<step*(t+1)-1;i++){
				runoff[t*step]=runoff[t*step]+runoff[i+1];	
			}
			sum[t]=runoff[t*step]/step;
		}

		return sum;

	}

}
