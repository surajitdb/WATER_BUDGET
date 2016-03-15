package waterBudget;

public class SimpleETModelFactory {

	public static EvapotranspirationModel createModel(String type,double ET,double S_i, double s_max){
		EvapotranspirationModel model=null;
		if (type.equals("AET")){
			model=new AETmodel(ET,S_i,s_max);
		}else if (type.equals("ExternalValues")){
			model=new ExternalETValues(ET);
		}
			
		return model;
		
	}
}
