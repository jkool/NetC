package conv;

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

public class ODAPExtract3D {

	private String conn = "http://tds.hycom.org/thredds/dodsC/glb_analysis";
	private String outputPath = "C:/Temp/INDx_2005.nc";
	private String inLatName = "Latitude";
	private String inLonName = "Longitude";
	private String inTimeName = "MT";
	private String inVarName = "ssh";
	private String outLatName = inLatName;
	private String outLonName = inLonName;
	private String outTimeName = "Time";
	private String outVarName = "ssh";
	private double minlon = 24;//Mnl -5.5;//IND 24;//SEAX 90;//IND 35;//NZ 160;//30;//74.2;//30;
	private double minlat = -40;//Mnl 30;//IND -40; //SEAX -20;//IND -30;//NZ -50;//-27.5;
	private int mintime = 37762;//Mnl 38352;//37762;//37562;//39643;//-1;// 37550;//
	private double maxlon = 106;//Mnl 42.5;//IND 106; //SEAX 175;//IND 100;//NZ 185;//74;//100;//74;
	private double maxlat = 30.5;//Mnl 47.5;//IND 30.5; // SEAX 40;//IND 31;//NZ -30;//31;
	private int maxtime = mintime+365;//Mnl 38716;//mintime+365;//Integer.MAX_VALUE;//39843;//Integer.MAX_VALUE;// 37700;//
	private int doover = 0;
	private int minx, miny, maxx, maxy, mint, maxt,full_length;
	private boolean neglon = false;
	private boolean split;
	private double terma, termb;
	private float[] ltr, lnr;
	private double[] tr;
	private Array lna, lta, ta;
	private ArrayList<Dimension> dims;
	NetcdfFileWriteable outfile;
	NetcdfDataset nds;

	public static void main(String[] args) {

		// If there are no arguments provided (e.g. file), try obtaining
		// from the keyboard.

		if (args == null) {

		}

		ODAPExtract3D to = new ODAPExtract3D();

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

		nds = NetcdfDataset.openDataset(conn);
		System.out.println(nds.toString());

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

	private static void getTextInput() {
		System.out.println("Implement me!");
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

		// Set up empty arrays for ranges and dimensions.

		dims = new ArrayList<Dimension>();

		// Create the output file

		outfile = NetcdfFileWriteable.createNew(outputPath, false);

		// Construct the data set dimensions - Time, Depth, Latitude and
		// Longitude (in order)

		Dimension timeDim = outfile.addDimension(outTimeName, tr.length);
		Dimension latDim = outfile.addDimension(outLatName, ltr.length);
		Dimension lonDim = outfile.addDimension(outLonName, lnr.length);

		// Add to a list - this becomes the co-ordinate system for the output
		// variable

		dims.add(timeDim);
		dims.add(latDim);
		dims.add(lonDim);

		// Create variables in the output file

		outfile.addVariable(outTimeName, DataType.DOUBLE, dims.subList(0, 1));
		outfile.addVariable(outLatName, DataType.DOUBLE, dims.subList(1, 2));
		outfile.addVariable(outLonName, DataType.DOUBLE, dims.subList(2, 3));
		outfile.addVariable(outVarName, DataType.FLOAT, dims);

		// Add attribute information (cloned from source)

		cloneAttributes(ncd, inTimeName, outfile, outTimeName);
		cloneAttributes(ncd, inLatName, outfile, outLatName);
		cloneAttributes(ncd, inLonName, outfile, outLonName);
		cloneAttributes(ncd, inVarName, outfile, outVarName);

		// Finalizes the structure of the output file, making the changes real.

		outfile.create();

		// Write the static information for 1D axes.

		outfile.write(outTimeName, ta);
		outfile.write(outLatName, lta);
		outfile.write(outLonName, lna);

		// Read the parameter variable.

		Variable v = ncd.findVariable(inVarName);

		// Rather than reading a 4-D chunk (which tends to be too large), we
		// instead loop across the time axis

		//for (int k = 0; k < maxt - mint + 1; k++) {
		for (int k = doover; k < maxt - mint; k++) {
		
			int m = 0;
			try {

				System.out.println("Converting time step " + k + "...");
				
				if(split){
					
					Array va = v.read(new int[] { k + mint, miny, minx },
							new int[] { 1, maxy - miny + 1,
							full_length - minx});
					
					outfile.write(outVarName, new int[] { k, 0, 0 }, va);
					
					Array vb = v.read(new int[] { k + mint, miny, 0 },
							new int[] { 1, maxy-miny + 1,
							maxx});
					
					outfile.write(outVarName, new int[] { k, 0, full_length-minx}, vb);
					
					}
					else{
				Array va = v.read(new int[] { k + mint, miny, minx },
						new int[] { 1, maxy - miny + 1, maxx - minx + 1 });
				outfile.write(outVarName, new int[] { k, 0, 0 }, va);}
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
							ncd = NetcdfDataset.openDataset(conn);
							repeat = false;
							nds = ncd;
							m = 0;
						} else if (resp.equalsIgnoreCase("n")) {
							System.out.println("Exiting.");
							System.exit(-1);
						} else {
							System.out
									.println(resp
											+ " is not a valid option.  Please select y(es) or n(o). > ");
						}
					}
				}
				try {
					Thread.sleep(5000 * (1 + m));
				} catch (InterruptedException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				m++;
				k--;
			}
		}

		// Clean-up duty

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

		Variable time = ncd.findVariable(inTimeName);

		if (time == null) {
			throw new IllegalArgumentException("Variable: " + inTimeName
					+ " was not found.");
		}

		// As a cheap start, we take the middle vector for both latitude and
		// longitude. For right now we're ducking the curvilinear coordinate
		// issue, and assuming we're not near the poles. We implement as a big
		// block instead of get1DDoubleRange, otherwise reading the array takes
		// forever. This is a consequence of using x and y as dimensions as
		// opposed to lat and lon. Lat and lon are not vectors, but 4-D arrays.

		lna = lon.read(new int[] { Math.round((lon.getShape()[0] - 1) / 2), 0,
				0, 0 }, new int[] { 1, lon.getShape()[1], 0, 0 });
		lta = lat.read(new int[] { 0, Math.round((lat.getShape()[1] - 1) / 2),
				0, 0 }, new int[] { lat.getShape()[0], 1, 0, 0 });
		ta = time.read();

		// Convert all of our 1D information into Java arrays for searching

		lnr = (float[]) lna.copyTo1DJavaArray();
		ltr = (float[]) lta.copyTo1DJavaArray();
		tr = (double[]) ta.copyTo1DJavaArray();

		// If we're using negative lon values as input, convert them to positive
		// values.

		if (neglon == true) {
			if (minlon < 0) {
				minlon = (360 + minlon)%360;
			}
			if (maxlon < 0) {
				maxlon = (360 + maxlon)%360;
			}
		}
		
		terma = (360+lnr[0])%360;
		termb = (360+lnr[lnr.length-1])%360;
		full_length = lnr.length;
		
		if(minlon<terma && maxlon>termb){split = true;}

		double df1 = Double.MAX_VALUE;
		double df2 = Double.MAX_VALUE;

		// Here we use modulo 360 - it is hard-coded in for the time being,
		// otherwise decision-making would be a giant pain to code. 
		// We find the minimum difference running through the entire data set.

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

		// Time (t) is not circular for our purposes, and
		// so instead we use binary search.

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
		
		if(split){
			float[] hold = new float[lnr.length-minx+maxx+1];
			System.arraycopy(lnr, minx, hold, 0, lnr.length-minx);
			System.arraycopy(lnr, 0, hold, lnr.length-minx, maxx+1);
			System.out.println();
			lnr = hold;
		}
		else{
		if (maxx - minx != lnr.length - 1) {
			lnr = Arrays.copyOfRange(lnr, minx, maxx + 1);
		}}
		if (maxy - miny != ltr.length - 1) {
			ltr = Arrays.copyOfRange(ltr, miny, maxy + 1);
		}
		if (maxt - mint != tr.length - 1) {
			tr = Arrays.copyOfRange(tr, mint, maxt);
		}

		lna = Array.factory(float.class, new int[] { lnr.length }, lnr);
		lta = Array.factory(float.class, new int[] { ltr.length }, ltr);
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
}

