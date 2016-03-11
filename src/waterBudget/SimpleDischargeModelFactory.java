package waterBudget;

public class SimpleDischargeModelFactory {

	public static DischargeModel createModel(String type,double Qinput, double A, 
			double a, double S_i, double b){
		DischargeModel model=null;
		if (type.equals("NonLinearReservoir")){
			model=new NonLinearReservoir(a,S_i,b);
		}else if (type.equals("ExternalValues")){
			model=new ExternalDischargeValues(Qinput,A);
		}
			
		return model;
		
	}
}
