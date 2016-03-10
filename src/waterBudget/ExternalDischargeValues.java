package waterBudget;

public class ExternalDischargeValues implements DischargeModel{
	
	double Qinput;
	double A;
	double Q;
	
	public ExternalDischargeValues(double Qinput, double A){
		this.A=A;
		this.Qinput=Qinput;
	}

	public double dischargeValues() {
		Q=Qinput/A*3.6;
		return Q;
	}
	
	

}
