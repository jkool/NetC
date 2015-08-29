package conv;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

import ucar.ma2.Array;
import ucar.ma2.ArrayDouble;
import ucar.ma2.DataType;
import ucar.ma2.Index;
import ucar.ma2.InvalidRangeException;
import ucar.nc2.Attribute;
import ucar.nc2.Dimension;
import ucar.nc2.NetcdfFileWriter;
import ucar.nc2.Variable;
import ucar.nc2.dataset.NetcdfDataset;
import utilities.MatrixUtilities;

/**
 * Extracts information from HYCOM using a network connection.
 * 
 * @author Johnathan Kool
 */

public class MakeTest_4D {

	private String outputPath = "C:\\Temp\\Zeros_21_v.nc";
	private String outLatName = "Latitude";
	private String outLonName = "Longitude";
	private String outLayerName = "Depth";
	private String outTimeName = "Time";
	private String outVarName = "v";
	private float minLon = -10;
	private float minLat = -10;
	private float minZ = -1;
	private double minTime = 0;
	private float maxLon = 10;
	private float maxLat = 10;
	private float maxZ = 0;
	private double maxTime = 100;
	
	private int latdim = 21;
	private int londim = 21;
	private int zdim = 21;
	private int tdim = 101;

	private float[] ltr, lnr, dr;
	private double[] tr;
	private ArrayList<Dimension> dims;
	NetcdfDataset nds;

	public static void main(String[] args) {

		// If there are no arguments provided (e.g. file), try obtaining
		// from the keyboard.

		if (args != null) {

		}

		MakeTest_4D mt = new MakeTest_4D();

		try {
			mt.run();
		} catch (IOException e) {
			// Error reading from file
			e.printStackTrace();
			System.exit(-1);
		}

		System.out.println("Complete.");
		System.exit(0);
	}

	private void run() throws IOException {

		try {
			outputSetup();
		} catch (InvalidRangeException e) {
			// Get text input and re-do ranges.
			System.out.println("Improper ranges provided.");
			e.printStackTrace();
			nds.close();
			System.exit(-1);
		}
	}

	/**
	 * This method sets up the output file for writing the output file. The
	 * potential mix of attribute types makes generalization complicated.
	 * 
	 * @param ncf
	 *            - The initial NetCDF file, used as a template for creating the
	 *            output file.
	 */

	private void outputSetup() throws InvalidRangeException,
			IOException {

		lnr = MatrixUtilities.linspace(minLon,maxLon,londim);
		ltr = MatrixUtilities.linspace(minLat,maxLat,latdim);
		dr = MatrixUtilities.linspace(maxZ,minZ,zdim);
		tr = MatrixUtilities.linspace(minTime,maxTime,tdim);
		
		// Set up empty arrays for ranges and dimensions.

		dims = new ArrayList<Dimension>();

		// Create the output file
		
		NetcdfFileWriter writer = NetcdfFileWriter.createNew(NetcdfFileWriter.Version.netcdf4, outputPath, null);
		
		// Construct the data set dimensions - Time, Depth, Latitude and
		// Longitude (in order)

		Dimension timeDim = writer.addUnlimitedDimension(outTimeName);
		Dimension layerDim = writer.addDimension(null,outLayerName, dr.length);
		Dimension latDim = writer.addDimension(null,outLatName, ltr.length);
		Dimension lonDim = writer.addDimension(null,outLonName, lnr.length);

		// Add to a list - this becomes the coordinate system for the output
		// variable

		dims.add(timeDim);
		dims.add(layerDim);
		dims.add(latDim);
		dims.add(lonDim);

		// Create variables in the output file

		Variable time = writer.addVariable(null,outTimeName, DataType.DOUBLE, Arrays.asList(timeDim));
		Variable layer = writer.addVariable(null,outLayerName, DataType.DOUBLE, Arrays.asList(layerDim));
		Variable lat = writer.addVariable(null,outLatName, DataType.DOUBLE, Arrays.asList(latDim));
		Variable lon = writer.addVariable(null,outLonName, DataType.DOUBLE, Arrays.asList(lonDim));

		//outfile.setLargeFile(true);
		Variable var = writer.addVariable(null,outVarName, DataType.FLOAT, dims);

		writer.addVariableAttribute(lat, new Attribute("units", "degrees_north"));
		writer.addVariableAttribute(lon, new Attribute("units", "degrees_east"));
		writer.addVariableAttribute(layer, new Attribute("units", "meters"));
		writer.addVariableAttribute(time, new Attribute("units", "units"));
		
		// Finalizes the structure of the output file, making the changes real.

		writer.create();

		// Write the static information for 1D axes.

		writer.write(time, Array.factory(tr));
		writer.write(layer, Array.factory(dr));
		writer.write(lat, Array.factory(ltr));
		writer.write(lon, Array.factory(lnr));

		ArrayDouble A = new ArrayDouble.D4(tdim, zdim, latdim, londim);
		Index idx = A.getIndex();
		
		// Loop across the time axis
		
		for (int t = 0; t < tdim; t++) {

			//Date d = new Date(System.currentTimeMillis());
			//System.out.println("Converting time step " + t + "..." + "("
			//		+ d.toString() + ")");
			
			//double [] val;
			//float ct = 0;
			
			for(int k = 0; k < zdim; k++){
				for(int i = 0; i < latdim; i++){
					for(int j = 0; j < londim; j++){
						A.set(idx.set(t,k,i,j), function(t,k,i,j));
					}
				}
			}
		}

		writer.write(var, new int[4], A);

		// Clean-up duty

		writer.flush();
		writer.close();
	}
	
	private double function(int t, int k, int i, int j){
		//if((20-j)<k){return Double.NaN;}return 1;
		return 0;
	}
	
	// Getters and setters

	public String getOutputPath() {
		return outputPath;
	}

	public void setOutputPath(String outputPath) {
		this.outputPath = outputPath;
	}

	public String getOutVarName() {
		return outVarName;
	}

	public void setOutVarName(String outVarName) {
		this.outVarName = outVarName;
	}

	public float getMinLon() {
		return minLon;
	}

	public void setMinLon(float minLon) {
		this.minLon = minLon;
	}

	public float getMinLat() {
		return minLat;
	}

	public void setMinLat(float minLat) {
		this.minLat = minLat;
	}

	public float getMinZ() {
		return minZ;
	}

	public void setMinZ(float minZ) {
		this.minZ = minZ;
	}

	public double getMinTime() {
		return minTime;
	}

	public void setMinTime(float minTime) {
		this.minTime = minTime;
	}

	public float getMaxLon() {
		return maxLon;
	}

	public void setMaxLon(float maxLon) {
		this.maxLon = maxLon;
	}

	public float getMaxLat() {
		return maxLat;
	}

	public void setMaxLat(float maxLat) {
		this.maxLat = maxLat;
	}

	public float getMaxZ() {
		return maxZ;
	}

	public void setMaxZ(float maxZ) {
		this.maxZ = maxZ;
	}

	public double getMaxTime() {
		return maxTime;
	}

	public void setMaxTime(float maxTime) {
		this.maxTime = maxTime;
	}
}
