package conv;

import ucar.nc2.*;
import ucar.ma2.*;

import java.io.*;
import java.util.*;

public class SquareCell2D {

	NetcdfFileWriteable ncw;
	NetcdfFile ncf;
	final float NODATA = 1E34f;
	String inLatName = "Latitude";
	String inLonName = "Longitude";
	String inDepthName = "depth";
	String inTimeName = "Time";
	String inVarName = "ssh";
	String outLatName = inLatName;
	String outLonName = inLonName;
	String outDepthName = inDepthName;
	String outTimeName = inTimeName;
	String outVarName = inVarName;
	String outVarDesc = "SSH";
	static String inFile = "C:/Temp/INDx_2005.nc";
	String outFile = "C:/Temp/INDx_2005_ssh_sq.nc";
	float lb_valid = -500f;
	float ub_valid = 500f;
	float minlat, maxlat, minlon, maxlon;

	public static void main(String[] args) {

		SquareCell2D sq = new SquareCell2D();
		try {
			sq.square(inFile);
		} catch (InvalidRangeException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.out.println("Complete.");
		System.exit(0);
	}

	private void square(String filename) throws InvalidRangeException,
			IOException {

		int nlat, nlon;

		ncf = NetcdfFile.open(filename);

		Variable lat = ncf.findVariable(inLatName);
		Variable lon = ncf.findVariable(inLonName);
		Variable time = ncf.findVariable(inTimeName);
		Variable var = ncf.findVariable(inVarName);

		int tlen = time.getShape()[0];

		Array latarr = lat.read();
		Array lonarr = lon.read();

		double[] platarr = (double[]) latarr.copyTo1DJavaArray();
		double[] plonarr = (double[]) lonarr.copyTo1DJavaArray();

		Index latindex = latarr.getIndex();
		Index lonindex = lonarr.getIndex();

		nlat = latarr.getShape()[0];
		nlon = lonarr.getShape()[0];

		minlat = latarr.getFloat(latindex.set(0));
		maxlat = latarr.getFloat(latindex.set(nlat - 1));
		float rangelat = maxlat - minlat;

		minlon = lonarr.getFloat(lonindex.set(0));
		maxlon = lonarr.getFloat(lonindex.set(nlon - 1));
		float rangelon = maxlon - minlon;

		// Keeping the same number of cells, what's the mean distance?

		float latdist = (rangelat) / (float) nlat;
		float londist = (rangelon) / (float) nlon;

		float cellsize = (latdist + londist) / 2f;

		int new_nlat = (int) Math.floor(rangelat / cellsize);
		int new_nlon = (int) Math.floor(rangelon / cellsize);

		// Create new arrays for use in dimensions

		float[] lataxis = new float[new_nlat];
		float[] lonaxis = new float[new_nlon];

		float latc = minlat;
		float lonc = minlon;

		for (int j = 0; j < new_nlat; j++) {
			lataxis[j] = latc;
			latc += cellsize;
		}

		for (int i = 0; i < new_nlon; i++) {
			lonaxis[i] = lonc;
			lonc += cellsize;
		}

		ncw = NetcdfFileWriteable.createNew(outFile, false);

		List<Dimension> dims = new ArrayList<Dimension>();

		dims.add(ncw.addDimension(outTimeName, tlen));
		dims.add(ncw.addDimension(outLatName, new_nlat));
		dims.add(ncw.addDimension(outLonName, new_nlon));

		ncw.addVariable(outTimeName, DataType.FLOAT, dims.subList(0, 1));

		ncw.addVariable(outLatName, DataType.FLOAT, dims.subList(1, 2));

		ncw.addVariable(outLonName, DataType.FLOAT, dims.subList(2, 3));

		ncw.addVariable(outVarName, DataType.FLOAT, dims);

		writeAttributeInformation();

		ncw.create();

		ncw.write(outTimeName, time.read());
		ncw.write(outLatName,
				Array.factory(float.class, new int[] { new_nlat }, lataxis));
		ncw.write(outLonName,
				Array.factory(float.class, new int[] { new_nlon }, lonaxis));

		float ilat, ilon;
		float ul, ll, ur, lr;
		double x, y, bs;
		int llat, llon, ulat, ulon;
		Array vararr;
		Array iarr;
		Index idx, odx;

		for (int t = 0; t < tlen; t++) {

			if (t % 10 == 0) {
				System.out.println("Processing step " + t + ". "
						+ (int) ((double) t / (double) tlen * 100)
						+ "% complete.");
			}

			vararr = var.read(new int[] { t, 0, 0 },
					new int[] { 1, nlat, nlon });
			vararr = vararr.reduce();
			idx = vararr.getIndex();
			iarr = Array.factory(float.class,
					new int[] { 1, new_nlat, new_nlon });
			odx = iarr.getIndex();

			for (int j = 0; j < new_nlat; j++) {
				for (int i = 0; i < new_nlon; i++) {

					// First we need the latitude and longitude of the point to
					// be interpolated.

					ilat = lataxis[j];
					ilon = lonaxis[i];

					llat = Arrays.binarySearch(platarr, ilat);
					if (llat < 0) {
						llat = -(llat + 2);
					}
					ulat = llat + 1;
					llon = Arrays.binarySearch(plonarr, ilon);
					if (llon < 0) {
						llon = -(llon + 2);
					}
					ulon = llon + 1;

					ul = vararr.getFloat(idx.set(ulat, llon));
					ur = vararr.getFloat(idx.set(ulat, ulon));
					ll = vararr.getFloat(idx.set(llat, llon));
					lr = vararr.getFloat(idx.set(llat, ulon));

					if (ul == NODATA || ur == NODATA || ll == NODATA
							|| lr == NODATA) {
						iarr.setFloat(odx.set(0, j, i), NODATA);

						// Bilinear interpolation

					} else {
						x = (ilon - plonarr[llon])
								/ (plonarr[ulon] - plonarr[llon]);
						y = (ilat - platarr[llat])
								/ (platarr[ulat] - platarr[llat]);
						bs = (ll * (1 - x) * (1 - y)) + (lr * x * (1 - y))
								+ (ul * (1 - x) * y) + (ur * x * y);
						iarr.setDouble(odx.set(0, j, i), bs);
					}
				}
			}

			ncw.write(outVarName, new int[] { t, 0, 0 }, iarr);

		}

		ncw.close();
	}

	/**
	 * Writes out the attribution information to the output file
	 */

	private void writeAttributeInformation() {

		cloneAttributes(ncf, inTimeName, ncw, outTimeName);
		ncw.addVariableAttribute(outLonName, "units", "degrees_east");
		ncw.addVariableAttribute(outLonName, "standard_name", outLonName);
		ncw.addVariableAttribute(outLonName, "axis", "x");
		ncw.addVariableAttribute(
				outLonName,
				"valid_range",
				Array.factory(float.class, new int[] { 2 }, new float[] {
						minlon, maxlon }));

		ncw.addVariableAttribute(outLatName, "units", "degrees_north");
		ncw.addVariableAttribute(outLatName, "standard_name", outLatName);
		ncw.addVariableAttribute(outLatName, "axis", "y");
		ncw.addVariableAttribute(
				outLatName,
				"valid_range",
				Array.factory(float.class, new int[] { 2 }, new float[] {
						minlat, maxlat }));
		ncw.addVariableAttribute(outVarName, "units", "ms-1");
		ncw.addVariableAttribute(outVarName, "missing_value", 1E34);
		ncw.addVariableAttribute(outVarName, "standard_name", outVarDesc);
		ncw.addVariableAttribute(outVarName, "long_name", outVarDesc);
		ncw.addVariableAttribute(
				outVarName,
				"valid_range",
				Array.factory(float.class, new int[] { 2 }, new float[] {
						lb_valid, ub_valid }));
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

	private void cloneAttributes(NetcdfFile infile, String inVarName,
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
