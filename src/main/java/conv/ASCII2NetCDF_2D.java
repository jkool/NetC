package conv;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.StringTokenizer;

import ucar.ma2.Array;
import ucar.ma2.ArrayFloat;
import ucar.ma2.DataType;
import ucar.ma2.Index;
import ucar.ma2.InvalidRangeException;
import ucar.nc2.Attribute;
import ucar.nc2.Dimension;
import ucar.nc2.NetcdfFileWriter;
import ucar.nc2.Variable;

public class ASCII2NetCDF_2D {

	private String variableName = "bathymetry";
	private String variableUnits = "meters";
	private String latName = "Latitude";
	private String lonName = "Longitude";
	private String latUnits = "degrees_north";
	private String lonUnits = "degrees_east";

	private int ncols;
	private int nrows;
	private float xllcorner;
	private float yllcorner;
	private float cellsize;
	private float NODATA = -9999;

	public void convert(String inputFile, String outputFile) {
		BufferedReader br = null;
		NetcdfFileWriter writer = null;

		try {

			// Read the header information

			br = new BufferedReader(new FileReader(inputFile));
			StringTokenizer stk;
			String ln = br.readLine();
			stk = new StringTokenizer(ln);
			stk.nextToken();
			ncols = Integer.parseInt(stk.nextToken());
			ln = br.readLine();
			stk = new StringTokenizer(ln);
			stk.nextToken();
			nrows = Integer.parseInt(stk.nextToken());
			ln = br.readLine();
			stk = new StringTokenizer(ln);
			stk.nextToken();
			xllcorner = Float.parseFloat(stk.nextToken());
			ln = br.readLine();
			stk = new StringTokenizer(ln);
			stk.nextToken();
			yllcorner = Float.parseFloat(stk.nextToken());
			ln = br.readLine();
			stk = new StringTokenizer(ln);
			stk.nextToken();
			cellsize = Float.parseFloat(stk.nextToken());
			ln = br.readLine();
			stk = new StringTokenizer(ln);
			stk.nextToken();
			NODATA = Integer.parseInt(stk.nextToken());

			// Generate the file

			writer = NetcdfFileWriter.createNew(NetcdfFileWriter.Version.netcdf4, outputFile, null);

			// Set up the dimensions

			float[] lats = new float[nrows];
			float[] lons = new float[ncols];

			for (int i = 0; i < nrows; i++) {
				lats[i] = yllcorner + (i * cellsize) + cellsize/2;
			}

			for (int i = 0; i < ncols; i++) {
				lons[i] = xllcorner + (i * cellsize) + cellsize/2;
			}

			Dimension latDim = writer.addDimension(null,latName, nrows);
			Dimension lonDim = writer.addDimension(null,lonName, ncols);

			ArrayList<Dimension> dims = new ArrayList<Dimension>();
			dims.add(latDim);
			dims.add(lonDim);

			// Add the Dimensions as variables

			Variable lat = writer.addVariable(null,latName, DataType.FLOAT,Arrays.asList(latDim));
			writer.addVariableAttribute(lat, new Attribute("units", latUnits));
			Variable lon = writer.addVariable(null, lonName, DataType.FLOAT, Arrays.asList(lonDim));
			writer.addVariableAttribute(lon, new Attribute("units", lonUnits));
			Variable var = writer.addVariable(null,variableName, DataType.FLOAT, dims);
			writer.addVariableAttribute(var, new Attribute("units", variableUnits));

			writer.create();

			Array lata = Array.factory(DataType.FLOAT,
					new int[] { lats.length }, lats);
			//lata.flip(0);
			Array lona = Array.factory(DataType.FLOAT,
					new int[] { lons.length }, lons);

			writer.write(lat, lata);
			writer.write(lon, lona);

			Array vardata = new ArrayFloat.D2(nrows, ncols);
			Index idx = vardata.getIndex();

			for (int i = nrows-1; i >=0; i--) {
				ln = br.readLine();
				stk = new StringTokenizer(ln);
				for (int j = 0; j < ncols; j++) {
					String s = stk.nextToken();
					float f = Float.parseFloat(s);
					vardata.setFloat(idx.set(i, j), f == NODATA ? Float.NaN : f);
				}
				System.out.format("Finished row %d of %d (%.2f%%)\n",nrows-i,nrows,100-(100*(float)i/nrows));
			}

			//vardata.flip(1);
			writer.write(var, vardata);
			writer.flush();
			writer.close();

		} catch (NumberFormatException e) {
			e.printStackTrace();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} catch (InvalidRangeException e) {
			e.printStackTrace();
		} finally {
			if (writer != null) {
				try {
					writer.close();
				} catch (IOException e) {
				}
			}
			if (br != null) {
				try {
					br.close();
				} catch (IOException e) {
				}
			}
		}
	}
	
	public static void main(String[] args){
		ASCII2NetCDF_2D a2 = new ASCII2NetCDF_2D();
		a2.convert("C:/Temp/lite.asc", "C:/Temp/aus_bath_lite2.nc");
	}
}
