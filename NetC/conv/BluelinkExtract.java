package conv;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;

import ucar.ma2.*;
import ucar.nc2.*;
import ucar.nc2.dataset.NetcdfDataset;

/**
 * Extracts information from HYCOM using a network connection.
 * 
 * @author Johnathan Kool
 */

public class BluelinkExtract {

	private String connpath = "http://www.marine.csiro.au/dods/nph-dods/dods-data/bl/BRAN2.1/u/";
	private String outputPath = "C:\\Temp\\BL-v.nc";
	private String inLatName = "yu_ocean";
	private String inLonName = "xu_ocean";
	private String inLayerName = "zt_ocean";
	private String inTimeName = "Time";
	private String inVarName = "u";
	private String outLatName = inLatName;
	private String outLonName = inLonName;
	private String outLayerName = "depth";
	private String outTimeName = "Time";
	private String outVarName = "u";
	private int minlon = 130;// 90
	private int minlat = -30;// -20
	private int mindpth = 0;
	private int mintime = -1;// 38250;
	private int maxlon = 175;// 170
	private int maxlat = 20;// 40
	private int maxdpth = 150;
	private int maxtime = 99999;// 38350;

	private int minx, miny, maxx, maxy, minz, maxz, mint, maxt;
	private boolean neglon = false;
	private float[] ltr, lnr, dr;
	private double[] tr;
	private Array lna, lta, da, ta;
	private ArrayList<Dimension> dims;
	NetcdfFileWriteable outfile;
	NetcdfDataset nds;

	public static void main(String[] args) {

		// If there are no arguments provided (e.g. file), try obtaining
		// from the keyboard.

		if (args == null) {

		}

		BluelinkExtract to = new BluelinkExtract();

		try {
			to.run();
		} catch (IOException e) {
			// Error reading from file
			e.printStackTrace();
			System.exit(-1);
		}

		System.out.println("Complete.");
		System.exit(0);

	}

	private void run() throws IOException {

		// Open the data set

		try {
			getBndBox(nds);
			outputSetup(nds);
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

	private void outputSetup(NetcdfDataset ncd) throws InvalidRangeException,
			IOException {

		FileReader fr = new FileReader("C://Temp//blink.txt");
		BufferedReader br = new BufferedReader(fr);
		
		ArrayList<String> slices = new ArrayList<String>();
		
		String ln = br.readLine();
		while(ln != null){
			slices.add(ln);
			ln = br.readLine();
		}
		
		// Set up empty arrays for ranges and dimensions.

		dims = new ArrayList<Dimension>();

		// Create the output file

		outfile = NetcdfFileWriteable.createNew(outputPath, false);

		// Construct the data set dimensions - Time, Depth, Latitude and
		// Longitude (in order)

		Dimension timeDim = outfile.addDimension(outTimeName, slices.size());
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

		outfile.addVariable(outTimeName, DataType.DOUBLE, dims.subList(0, 1));
		outfile.addVariable(outLayerName, DataType.DOUBLE, dims.subList(1, 2));
		outfile.addVariable(outLatName, DataType.DOUBLE, dims.subList(2, 3));
		outfile.addVariable(outLonName, DataType.DOUBLE, dims.subList(3, 4));

		// outfile.setLargeFile(true);
		outfile.addVariable(outVarName, DataType.FLOAT, dims);

		// Add attribute information (cloned from source)

		cloneAttributes(ncd, inTimeName, outfile, outTimeName);
		cloneAttributes(ncd, inLayerName, outfile, outLayerName);
		cloneAttributes(ncd, inLatName, outfile, outLatName);
		cloneAttributes(ncd, inLonName, outfile, outLonName);
		cloneAttributes(ncd, inVarName, outfile, outVarName);

		// Finalizes the structure of the output file, making the changes real.

		outfile.create();

		// Write the static information for 1D axes.

		outfile.write(outTimeName, ta);
		outfile.write(outLayerName, da);
		outfile.write(outLatName, lta);
		outfile.write(outLonName, lna);

		// Read the parameter variable.

		int m = 0;

		for (int k = 0; k < slices.size(); k++) {

		ncd = NetcdfDataset.openDataset(connpath + slices.get(k));
		Variable v = ncd.findVariable(inVarName);

		// Rather than reading a 4-D chunk (which tends to be too large), we
		// instead loop across the time axis

			System.out.println("Converting time step " + k + "...");
			try {
				Array va = v.read(new int[] { k + mint, minz, miny, minx },
						new int[] { 1, maxz - minz + 1, maxy - miny + 1,
								maxx - minx + 1 });
				outfile.write(outVarName, new int[] { k, 0, 0, 0 }, va);
			} catch (Exception de) {
				System.out.println("WARNING:  OPENDAP Error - retry");
				if (m > 5) {
					Scanner input = new Scanner(System.in);
					System.out
							.println("Exceeded the set number of retries.  Continue?  (y/n) > ");
					boolean repeat = true;
					while (repeat) {
						String resp = input.next();
						if (resp.equalsIgnoreCase("y")) {
							ncd.close();
							ncd = NetcdfDataset.openDataset(connpath + slices.get(k));
							m = 0;
							repeat = false;
							nds = ncd;
						} else if (resp.equalsIgnoreCase("n")) {
							System.out.println("Exiting");
							System.exit(-1);
						}
						else{
							System.out.println(resp + " is not a valid option.  Please select y(es) or n(o). > ");
						}
					}
					input.close();
				}
				try {
					Thread.sleep(5000);
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				m++;
				k--;
			}
		}

		// Clean-up duty

		br.close();
		outfile.flush();
		outfile.finish();
		outfile.close();
	}

	/**
	 * Clips out the region of interest (spatial and temporal)
	 * 
	 * @param ncd
	 * @throws InvalidRangeException
	 * @throws IOException
	 */

	private void getBndBox(NetcdfDataset ncd) throws InvalidRangeException,
			IOException {

		// Read in all the variables.

		Variable lon = ncd.findVariable(inLonName);

		if (lon == null) {
			throw new IllegalArgumentException("Variable: " + inLonName
					+ " was not found.");
		}

		Variable lat = ncd.findVariable(inLatName);

		if (lat == null) {
			throw new IllegalArgumentException("Variable: " + inLatName
					+ " was not found.");
		}

		Variable layer = ncd.findVariable(inLayerName);

		if (layer == null) {
			throw new IllegalArgumentException("Variable: " + inLayerName
					+ " was not found.");
		}

		Variable time = ncd.findVariable(inTimeName);

		if (time == null) {
			throw new IllegalArgumentException("Variable: " + inTimeName
					+ " was not found.");
		}

		// There is a possibility that the reference scheme may be shifted.
		// Also, since this is circular let's do a brute force search and
		// look for the minimum distance between values.

		// First, we need to retrieve index vectors for lon and lat. This is
		// done by reading from lon and lat, the first term is the origin, and
		// the second is the shape we wish to retrieve.
		// As a cheap start, we take the middle vector for both latitude and
		// longitude. For right now we're ducking the curvilinear coordinate
		// issue, and assuming we're not near the poles. We implement as a big
		// block instead of get1DDoubleRange, otherwise reading the array takes
		// forever. This is a consequence of using x and y as dimensions as
		// opposed
		// to lat and lon. Lat and lon are no longer vectors, but 4-D arrays.

		lna = lon.read();
		// lna = lon.read(new int[] { Math.round((lon.getShape()[0] - 1) / 2),
		// 0}, new int[] { 1, lon.getShape()[1]});
		lta = lat.read();
		// lta = lat.read(new int[] { 0, Math.round((lat.getShape()[1] - 1) /
		// 2)}, new int[] { lat.getShape()[0], 1});
		da = layer.read();
		ta = time.read();

		// Convert all of our 1D information into Java arrays for searching

		lnr = (float[]) lna.copyTo1DJavaArray();
		ltr = (float[]) lta.copyTo1DJavaArray();
		dr = (float[]) da.copyTo1DJavaArray();
		tr = (double[]) ta.copyTo1DJavaArray();

		// If we're using negative lon values as input, convert them to positive
		// values.

		if (neglon == true) {
			if (minlon < 0) {
				minlon = 360 + minlon;
			}
			if (maxlon < 0) {
				maxlon = 360 + maxlon;
			}
		}

		double df1 = Double.MAX_VALUE;
		double df2 = Double.MAX_VALUE;

		// Here we use modulo 360 - it is hard-coded in for the time being,
		// otherwise
		// decision-making would be a giant pain to code. We find the minimum
		// difference
		// running through the entire data set.

		for (int i = 0; i < lnr.length; i++) {
			lnr[i] = lnr[i] % 360;
			if (Math.abs(lnr[i] - minlon) < Math.abs(df1)) {
				df1 = lnr[i] - minlon;
				minx = i;
			}
			if (Math.abs(lnr[i] - maxlon) < Math.abs(df2)) {
				df2 = lnr[i] - maxlon;
				maxx = i;
			}
		}

		df1 = Double.MAX_VALUE;
		df2 = Double.MAX_VALUE;

		for (int i = 0; i < ltr.length; i++) {
			if (Math.abs(ltr[i] - minlat) < Math.abs(df1)) {
				df1 = ltr[i] - minlat;
				miny = i;
			}
			if (Math.abs(ltr[i] - maxlat) < Math.abs(df2)) {
				df2 = ltr[i] - maxlat;
				maxy = i;
			}
		}

		// Depth (z) and Time (t) are not circular for our purposes, and
		// so instead we use binary search.

		minz = Arrays.binarySearch(dr, mindpth);
		if (minz < 0) {
			minz = -(minz + 1);
		}

		maxz = Arrays.binarySearch(dr, maxdpth);
		if (maxz < 0) {
			maxz = -(maxz + 1);
		}

		mint = Arrays.binarySearch(tr, mintime);
		if (mint < 0) {
			mint = -(mint + 1);
		}

		maxt = Arrays.binarySearch(tr, maxtime);
		if (maxt < 0) {
			if (maxt == -1) {
				maxt = tr.length - 1;
			} else {

				maxt = -(maxt + 1);
			}
		}

		// Since we're here and have all the variables we need, we might as well
		// clip now

		if (maxx - minx != lnr.length - 1) {
			lnr = Arrays.copyOfRange(lnr, minx, maxx + 1);
		}
		if (maxy - miny != ltr.length - 1) {
			ltr = Arrays.copyOfRange(ltr, miny, maxy + 1);
		}
		if (maxz - minz != dr.length - 1) {
			dr = Arrays.copyOfRange(dr, minz, maxz + 1);
		}
		if (maxt - mint != tr.length - 1) {
			tr = Arrays.copyOfRange(tr, mint, maxt);
		}

		lna = Array.factory(float.class, new int[] { lnr.length }, lnr);
		lta = Array.factory(float.class, new int[] { ltr.length }, ltr);
		da = Array.factory(float.class, new int[] { dr.length }, dr);
		ta = Array.factory(double.class, new int[] { tr.length }, tr);

		for (int i = 0; i < tr.length - 1; i++) {
			if (tr[i + 1] - tr[i] != 1) {
				System.out.println(i + ": " + tr[i] + ", " + tr[i + 1]);
			}
		}
	}

	/**
	 * Clones attribute information from source variables to destination
	 * variables.
	 * 
	 * @param infile
	 *            - The input file
	 * @param inVarName
	 *            - The name of the input variable
	 * @param outfile
	 *            - The output file
	 * @param outVarName
	 *            - The name of the output variable
	 */

	private void cloneAttributes(NetcdfDataset infile, String inVarName,
			NetcdfFileWriteable outfile, String outVarName) {

		// Find the variable

		Variable vi = infile.findVariable(inVarName);

		// Grab all of its attributes - unchecked, but should be OK.

		List<Attribute> l = vi.getAttributes();
		for (Attribute a : l) {

			outfile.addVariableAttribute(outVarName, a);
		}
	}
	
	public String getOutputPath() {
		return outputPath;
	}

	public void setOutputPath(String outputPath) {
		this.outputPath = outputPath;
	}

	public String getInVarName() {
		return inVarName;
	}

	public void setInVarName(String inVarName) {
		this.inVarName = inVarName;
	}

	public String getOutVarName() {
		return outVarName;
	}

	public void setOutVarName(String outVarName) {
		this.outVarName = outVarName;
	}

	public int getMinlon() {
		return minlon;
	}

	public void setMinlon(int minlon) {
		this.minlon = minlon;
	}

	public int getMinlat() {
		return minlat;
	}

	public void setMinlat(int minlat) {
		this.minlat = minlat;
	}

	public int getMindpth() {
		return mindpth;
	}

	public void setMindpth(int mindpth) {
		this.mindpth = mindpth;
	}

	public int getMintime() {
		return mintime;
	}

	public void setMintime(int mintime) {
		this.mintime = mintime;
	}

	public int getMaxlon() {
		return maxlon;
	}

	public void setMaxlon(int maxlon) {
		this.maxlon = maxlon;
	}

	public int getMaxlat() {
		return maxlat;
	}

	public void setMaxlat(int maxlat) {
		this.maxlat = maxlat;
	}

	public int getMaxdpth() {
		return maxdpth;
	}

	public void setMaxdpth(int maxdpth) {
		this.maxdpth = maxdpth;
	}

	public int getMaxtime() {
		return maxtime;
	}

	public void setMaxtime(int maxtime) {
		this.maxtime = maxtime;
	}
}
