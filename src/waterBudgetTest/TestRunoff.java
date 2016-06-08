package waterBudgetTest;


import java.net.URISyntaxException;
import java.util.HashMap;

import org.geotools.coverage.grid.GridCoverage2D;
import org.jgrasstools.gears.io.rasterreader.OmsRasterReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorWriter;
import org.jgrasstools.hortonmachine.utils.HMTestCase;

import runoff.WaterBudget;

public class TestRunoff extends HMTestCase{

	public void testLinear() throws Exception {

		String startDate = "1994-01-01 00:00";
		String endDate = "1994-01-02 00:00";
		int timeStepMinutes = 60;
		String fId = "ID";

		String inPathToPrec = "resources/Input/rainfall.csv";
		String pathToQ= "resources/Output/runoff/Q_runoff.csv";


		
		OmsTimeSeriesIteratorReader JReader = getTimeseriesReader(inPathToPrec, fId, startDate, endDate, timeStepMinutes);

		OmsTimeSeriesIteratorWriter writerQ = new OmsTimeSeriesIteratorWriter();

		
		writerQ.file = pathToQ;
		writerQ.tStart = startDate;
		writerQ.tTimestep = timeStepMinutes;
		writerQ.fileNovalue="-9999";
		
		

		
		WaterBudget waterBudget= new WaterBudget();
		
		OmsRasterReader Wsup = new OmsRasterReader();
		Wsup.file = "resources/Input/width_10.asc";
		Wsup.fileNovalue = -9999.0;
		Wsup.geodataNovalue = Double.NaN;
		Wsup.process();
		GridCoverage2D width_sup = Wsup.outRaster;
		
		
		OmsRasterReader topindex = new OmsRasterReader();
		topindex.file = "resources/Input/topIndex.asc";
		topindex.fileNovalue = -9999.0;
		topindex.geodataNovalue = Double.NaN;
		topindex.process();
		GridCoverage2D topIndex = topindex.outRaster;


		while( JReader.doProcess ) {
		
			waterBudget.solver_model="dp853";
			waterBudget.inRescaledDistance=width_sup;
			waterBudget.pCelerity=2;
			waterBudget.alpha=0.5;
			waterBudget.inTopindex=topIndex;
			waterBudget.pSat=2;
			waterBudget.inTimestep=timeStepMinutes;
			waterBudget.tStartDate=startDate;
			waterBudget.tEndDate=endDate;
			waterBudget.ID=209;
			
			JReader.nextRecord();
			
			HashMap<Integer, double[]> id2ValueMap = JReader.outData;
			waterBudget.inRainValues = id2ValueMap;

            waterBudget.process();
            
            HashMap<Integer, double[]> outHMDischarge = waterBudget.outHMDischarge;

			
			writerQ.inData = outHMDischarge;
			writerQ.writeNextLine();
			
			if (pathToQ != null) {
				writerQ.close();
			}
          
		}
		JReader.close();


	}


	private OmsTimeSeriesIteratorReader getTimeseriesReader( String inPath, String id, String startDate, String endDate,
			int timeStepMinutes ) throws URISyntaxException {
		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file = inPath;
		reader.idfield = "ID";
		reader.tStart = startDate;
		reader.tTimestep = timeStepMinutes;
		reader.tEnd = endDate;
		reader.fileNovalue = "-9999";
		reader.initProcess();
		return reader;
	}
}