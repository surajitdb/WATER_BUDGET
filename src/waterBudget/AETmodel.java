package waterBudget;

public class AETmodel implements EvapotranspirationModel{

	public double s_max;
	public double ET;
	public double S_i;
	public double AET;

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
