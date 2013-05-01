package conv;

import java.io.IOException;
import java.util.ArrayList;
//import java.util.Date;

import ucar.ma2.*;
import ucar.nc2.*;
import ucar.nc2.dataset.NetcdfDataset;
import utilities.MatrixUtilities;

/**
 * Extracts information from HYCOM using a network connection.
 * 
 * @author Johnathan Kool
 */

public class MakeTest_2D {

	private String outputPath = "C:\\Temp\\Floor_increasing_x_slope.nc";
	private String outLatName = "Latitude";
	private String outLonName = "Longitude";
	private String outVarName = "bathymetry";
	private float minLon = -10;
	private float minLat = -10;
	private float maxLon = 10;
	private float maxLat = 10;

	private int latdim = 21;
	private int londim = 21;

	private float[] ltr, lnr;
	private ArrayList<Dimension> dims;
	NetcdfFileWriteable outfile;
	NetcdfDataset nds;

	public static void main(String[] args) {

		// If there are no arguments provided (e.g. file), try obtaining
		// from the keyboard.

		if (args != null) {

		}

		MakeTest_2D mt = new MakeTest_2D();

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

	private void outputSetup() throws InvalidRangeException, IOException {

		lnr = MatrixUtilities.linspace(minLon, maxLon, londim);
		ltr = MatrixUtilities.linspace(minLat, maxLat, latdim);

		// Set up empty arrays for ranges and dimensions.

		dims = new ArrayList<Dimension>();

		// Create the output file

		outfile = NetcdfFileWriteable.createNew(outputPath, false);

		// Construct the data set dimensions - Time, Depth, Latitude and
		// Longitude (in order)

		Dimension latDim = outfile.addDimension(outLatName, ltr.length);
		Dimension lonDim = outfile.addDimension(outLonName, lnr.length);

		// Add to a list - this becomes the coordinate system for the output
		// variable

		dims.add(latDim);
		dims.add(lonDim);

		// Create variables in the output file

		outfile.addVariable(outLatName, DataType.DOUBLE,
				new Dimension[] { latDim });
		outfile.addVariable(outLonName, DataType.DOUBLE,
				new Dimension[] { lonDim });

		// outfile.setLargeFile(true);
		outfile.addVariable(outVarName, DataType.FLOAT, dims);

		outfile.addVariableAttribute(outLatName, "units", "degrees_north");
		outfile.addVariableAttribute(outLonName, "units", "degrees_east");

		// Finalizes the structure of the output file, making the changes real.

		outfile.create();

		// Write the static information for 1D axes.

		outfile.write(outLatName, Array.factory(ltr));
		outfile.write(outLonName, Array.factory(lnr));

		ArrayDouble A = new ArrayDouble.D2(latdim, londim);
		Index idx = A.getIndex();

		//double[] val;
		//float ct = 0;
		for (int i = 0; i < latdim; i++) {
			for (int j = 0; j < londim; j++) {
				// float scale = .00001f/1.11f;
				A.set(idx.set(i, j), function(i,j));
				// A.set(idx.set(i,j), ltr[i]);
				// A.set(idx.set(t,k,i,j), ct*1E-3);
				// val = Utils.lonlat2ceqd(new double[]{ltr[i],lnr[j]});
				// A.set(idx.set(i,j), (float)val[0]);
				// A.set(idx.set(i,j), ct);
				//ct++;
			}
		}

		outfile.write(outVarName, new int[2], A);

		// Clean-up duty

		outfile.flush();
		outfile.finish();
		outfile.close();
	}

	// Getters and setters
	
	public double function(int i, int j){
		//return -1+((double)j/21d);
		return Math.max((i-10.5)*(i-10.5), (j-10.5)*(j-10.5));
	}
	
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
}
