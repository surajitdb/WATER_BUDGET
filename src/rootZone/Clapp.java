package rootZone;

public class Clapp implements DischargeModel{
	
	double a;
	double S_i;
	double b;
	double Q;
	
	public Clapp(double a,double S_i, double b){
		this.a=a;
		this.S_i=S_i;
		this.b=b;		
	}

	public double dischargeValues() {
		Q=a*Math.pow(S_i, 2*b+3);
		return Q;
	}

}
