package conv;

import java.io.IOException;
import java.util.ArrayList;
//import java.util.Date;

import ucar.ma2.*;
import ucar.nc2.*;
import ucar.nc2.dataset.NetcdfDataset;
import utilities.MatrixUtilities;
//import utilities.Utils;

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
	NetcdfFileWriteable outfile;
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
		
		outfile = NetcdfFileWriteable.createNew(outputPath, false);
		
		// Construct the data set dimensions - Time, Depth, Latitude and
		// Longitude (in order)

		Dimension timeDim = outfile.addUnlimitedDimension(outTimeName);
		Dimension layerDim = outfile.addDimension(outLayerName, dr.length);
		Dimension latDim = outfile.addDimension(outLatName, ltr.length);
		Dimension lonDim = outfile.addDimension(outLonName, lnr.length);

		// Add to a list - this becomes the coordinate system for the output
		// variable

		dims.add(timeDim);
		dims.add(layerDim);
		dims.add(latDim);
		dims.add(lonDim);

		// Create variables in the output file

		outfile.addVariable(outTimeName, DataType.DOUBLE, new Dimension[]{timeDim});
		outfile.addVariable(outLayerName, DataType.DOUBLE, new Dimension[]{layerDim});
		outfile.addVariable(outLatName, DataType.DOUBLE, new Dimension[]{latDim});
		outfile.addVariable(outLonName, DataType.DOUBLE, new Dimension[]{lonDim});

		//outfile.setLargeFile(true);
		outfile.addVariable(outVarName, DataType.FLOAT, dims);

		outfile.addVariableAttribute(outLatName, "units", "degrees_north");
		outfile.addVariableAttribute(outLonName, "units", "degrees_east");
		outfile.addVariableAttribute(outLayerName, "units", "meters");
		outfile.addVariableAttribute(outTimeName, "units", "units");
		
		// Finalizes the structure of the output file, making the changes real.

		outfile.create();

		// Write the static information for 1D axes.

		outfile.write(outTimeName, Array.factory(tr));
		outfile.write(outLayerName, Array.factory(dr));
		outfile.write(outLatName, Array.factory(ltr));
		outfile.write(outLonName, Array.factory(lnr));

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
						//val =  -((2*dr[k])-1);
						//float[] arr = new float[]{lnr[j],ltr[i],dr[k]};
						//float[] arr = new float[]{lnr[j],ltr[i]};
						//float[] tmp = MatrixUtilities.negative(MatrixUtilities.cross(arr, new float[]{1f,1f,1f}));
						//float [] n = MatrixUtilities.normalize(arr);
						//float[] gravity = MatrixUtilities.subtract(n,arr);
						//val = MatrixUtilities.add(tmp, gravity);
						//float scale = .00001f/1.11f;
						A.set(idx.set(t,k,i,j), function(t,k,i,j));
						//A.set(idx.set(t,k,i,j), val[0]);
						//A.set(idx.set(t,k,i,j), val[1]);
						//A.set(idx.set(t,k,i,j), dr[k]>=0.5?Float.NaN:dr[k]);
						//A.set(idx.set(t,k,i,j), dr[k]>=0.5?Float.NaN:-ltr[i]*scale);
						//A.set(idx.set(t,k,i,j), ct*1E-3);
						//val = Utils.lonlat2ceqd(new double[]{ltr[i],lnr[j]});
						//A.set(idx.set(t,k,i,j), lnr[j]);
						//A.set(idx.set(t,k,i,j), dr[k]>=0.5?Float.NaN:0);
						//ct++;
					}
				}
			}
		}

		outfile.write(outVarName, new int[4], A);

		// Clean-up duty

		outfile.flush();
		outfile.finish();
		outfile.close();
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
