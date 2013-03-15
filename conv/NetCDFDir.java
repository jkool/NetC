package conv;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

import ucar.ma2.Array;
import ucar.nc2.NetcdfFile;
import ucar.nc2.Variable;

public class NetCDFDir {

	private Map<Long, NetcdfFile> files = new TreeMap<Long, NetcdfFile>();
	private String tName = "MT";
	
	public NetCDFDir(String s){
		this(s,new FilenamePatternFilter(".*"));
	}

	public NetCDFDir(String s, FilenameFilter filter) {

		File[] fa = new File(s).listFiles(filter);

		try {
			for (File fil : fa) {

				NetcdfFile ncf = NetcdfFile.open(fil.getPath());
				Variable tVar = ncf.findVariable(tName);

				Array arr = tVar.read();

				// Convert into a java array

				double[] ja = (double[]) arr.copyTo1DJavaArray();

				// Determine the minimum and maximum times

				long[] minmax = new long[2];

				minmax[0] = TimeConvert.HYCOMToMillis((long) ja[0]);
				minmax[1] = TimeConvert.HYCOMToMillis((long) ja[ja.length - 1]);

				// Put into an index linking start time with the associated file

			    files.put(minmax[0], ncf);				
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public Map<Long, NetcdfFile> getFiles() {
		return files;
	}
	
	public Set<Long> getTimes(){
		return files.keySet();
	}

	public void setFiles(Map<Long, NetcdfFile> files) {
		this.files = files;
	}

	public String gettName() {
		return tName;
	}
	
	/**
	 * Sets the name of the Time variable
	 * 
	 * @param tName
	 */

	public void settName(String tName) {
		this.tName = tName;
	}
}
