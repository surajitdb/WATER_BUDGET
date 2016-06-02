package groundWater;

public class AETmodel implements ETModel{
	
	public static double s_max;
	public static double ET;
	public static double S_i;
	public static double AET;

	
	public AETmodel(double ET,double S_i, double s_max){
		this.s_max=s_max;
		this.S_i=S_i;
		this.ET=ET;
	}

	public double ETValues() {
		AET=ET*S_i/s_max;
		return AET;
	}



}
