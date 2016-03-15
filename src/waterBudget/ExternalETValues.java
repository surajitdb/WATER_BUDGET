package waterBudget;

public class ExternalETValues implements EvapotranspirationModel{
	
	double ETinput;
	double ET;
	
	public ExternalETValues(double ETinput){
		this.ETinput=ETinput;
	}


	public double ETValues() {
		ET=ETinput;
		return ET;
	}
	
	

}
