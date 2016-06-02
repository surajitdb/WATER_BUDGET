package rootZone;

public class SimpleETModelFactory {

	public static ETModel createModel(String type,double ET,double S_i, double s_max){
		ETModel model=null;
		if (type.equals("AET")){
			model=new AETmodel(ET,S_i,s_max);
		}else if (type.equals("ExternalValues")){
			model=new ExternalETValues(ET);
		}else if (type.equals("nullET")){
			model=new AETnull();
		}
			
		return model;
		
	}
}
