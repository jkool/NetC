package conv;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import ucar.ma2.Array;
import ucar.ma2.DataType;
import ucar.ma2.Index;
import ucar.ma2.InvalidRangeException;
import ucar.nc2.Attribute;
import ucar.nc2.Dimension;
import ucar.nc2.NetcdfFileWriter;
import ucar.nc2.Variable;
import ucar.nc2.dataset.NetcdfDataset;
import utilities.MatrixUtilities;

public class MakeW {

	// public static void main(String[] args) {
	// MakeW mw = new MakeW();
	// mw.go();
	// }
	private String inputUFile = "D:/HYCOM/Blocks/2008_12_30_to_2009_01_29_tmp_u.nc";
	private String inputVFile = "D:/HYCOM/Blocks/2008_12_30_to_2009_01_29_tmp_v.nc";
	private String outputWFile = "D:/HYCOM/Blocks/2008_12_30_to_2009_01_29_tmp_w.nc";

	private ArrayList<Dimension> dims;
	private NetcdfDataset uFile, vFile;
	NetcdfFileWriter writer;
	private String inLatName = "Latitude";
	private String inLonName = "Longitude";
	private String inLayerName = "Depth";
	private String inTimeName = "Time";
	private Variable var;
	private String inUName = "u";
	private String inVName = "v";
	private String outLatName = inLatName;
	private String outLonName = inLonName;
	private String outLayerName = "Depth";
	private String outTimeName = "Time";
	private String outVarName = "w";
	private boolean reproject = true;

	/**
	 * Calculates the change in w over a change in z
	 * 
	 * @param lats
	 *            - an Array of latitude values
	 * @param lons
	 *            - an Array of longitude values
	 * @param us
	 *            - an Array of u (east-west) velocity values
	 * @param vs
	 *            - an Array of v (north-south) velocity values
	 * @return
	 */

	public float calcdwdz(Array lats, Array lons, Array us, Array vs) {
		float du = dx(us);
		if (Float.isNaN(du)) {
			return du;
		}
		float dv = dy(vs);
		if (Float.isNaN(dv)) {
			return dv;
		}

		Array[] prj;

		if (reproject) {
			prj = prj2meters(lons, lats);
		} else {
			prj = new Array[] { lons, lats };
		}

		float dx = dx(prj[1]);
		float dy = dy(prj[0]);
		float dudx = du / dx;
		float dvdy = dv / dy;
		float out = -(dudx + dvdy);
		return out==-0.0f?0.0f:out;
	}

	/**
	 * Clones the attributes of the input file(s) to the output file
	 * 
	 * @param infile
	 * @param inVarName
	 * @param outfile
	 * @param outVarName
	 */

	private void cloneAttributes(NetcdfDataset infile, String inVarName,
			NetcdfFileWriter writer, Variable outVar) {

		// Find the variable

		Variable vi = infile.findVariable(inVarName);

		List<Attribute> l = vi.getAttributes();

		for (Attribute a : l) {
			writer.addVariableAttribute(outVar, a);
		}
	}

	/**
	 * Calculates change in the x direction
	 * 
	 * @param arr
	 *            - a 3x3 array of numbers
	 * @return
	 */

	public float dx(Array arr) {
		Index idx = Index.factory(arr.getShape());
		
		float e = arr.getFloat((int) arr.getSize()/2);
		if (Float.isNaN(e)) {
			return Float.NaN;
		}
		
		if(arr.getSize()==3){
			return ((arr.getFloat(2)-arr.getFloat(0))/2);
		}
		
		float a = arr.getFloat(idx.set(0, 0));
		float c = arr.getFloat(idx.set(0, 2));
		float d = arr.getFloat(idx.set(1, 0));
		float f = arr.getFloat(idx.set(1, 2));
		float g = arr.getFloat(idx.set(2, 0));
		float i = arr.getFloat(idx.set(2, 2));
		a = Float.isNaN(a) ? 0 : a;
		c = Float.isNaN(c) ? 0 : c;
		d = Float.isNaN(d) ? 0 : d;
		f = Float.isNaN(f) ? 0 : f;
		g = Float.isNaN(g) ? 0 : g;
		i = Float.isNaN(i) ? 0 : i;
		return ((c + 2 * f + i) - (a + 2 * d + g)) / 8;
	}

	/**
	 * Calculates change in the x direction
	 * 
	 * @param arr1
	 *            , arr2 - 3x3 arrays of numbers
	 * @return
	 */

	/**
	 * Calculates change in the y direction
	 * 
	 * @param arr
	 *            - a 3x3 array of numbers
	 * @return
	 */

	public float dy(Array arr) {
		Index idx = Index.factory(arr.getShape());
		
		float e = arr.getFloat((int) arr.getSize()/2);
		if (Float.isNaN(e)) {
			return Float.NaN;
		}
		
		if(arr.getSize()==3){
			return ((arr.getFloat(2)-arr.getFloat(0))/2);
		}
		
		float a = arr.getFloat(idx.set(0, 0));
		float b = arr.getFloat(idx.set(0, 1));
		float c = arr.getFloat(idx.set(0, 2));
		float g = arr.getFloat(idx.set(2, 0));
		float h = arr.getFloat(idx.set(2, 1));
		float i = arr.getFloat(idx.set(2, 2));
		a = Float.isNaN(a) ? 0 : a;
		b = Float.isNaN(b) ? 0 : b;
		c = Float.isNaN(c) ? 0 : c;
		g = Float.isNaN(g) ? 0 : g;
		h = Float.isNaN(h) ? 0 : h;
		i = Float.isNaN(i) ? 0 : i;
		return ((g + 2 * h + i) - (a + 2 * b + c)) / 8;
	}

	/**
	 * Creates the w file
	 * 
	 * @param uFile
	 *            - the input east-west velocity file
	 * @param vFile
	 *            - the input north-south velocity file
	 * @throws InvalidRangeException
	 * @throws IOException
	 */

	public NetcdfDataset generate(NetcdfDataset uFile, NetcdfDataset vFile)
			throws InvalidRangeException, IOException {

		// Set up empty arrays for ranges and dimensions.

		dims = new ArrayList<Dimension>();

		// Create the output file

		writer = NetcdfFileWriter.createNew(NetcdfFileWriter.Version.netcdf4, outputWFile, null);

		// Construct the data set dimensions - Time, Depth, Latitude and
		// Longitude (in order)

		int uTimeLength = uFile.findVariable(inTimeName).getShape()[0];
		int uLayerLength = uFile.findVariable(inLayerName).getShape()[0];
		int uLatLength = uFile.findVariable(inLatName).getShape()[0];
		int uLonLength = uFile.findVariable(inLonName).getShape()[0];

		int vTimeLength = vFile.findVariable(inTimeName).getShape()[0];
		int vLayerLength = vFile.findVariable(inLayerName).getShape()[0];
		int vLatLength = vFile.findVariable(inLatName).getShape()[0];
		int vLonLength = vFile.findVariable(inLonName).getShape()[0];

		Dimension timeDim = writer.addDimension(null, outTimeName,
				Math.min(uTimeLength, vTimeLength));
		Dimension layerDim = writer.addDimension(null, outLayerName,
				Math.min(uLayerLength, vLayerLength));
		Dimension latDim = writer.addDimension(null, outLatName,
				Math.min(uLatLength, vLatLength));
		Dimension lonDim = writer.addDimension(null, outLonName,
				Math.min(uLonLength, vLonLength));

		// Add to a list - this becomes the coordinate system for the output
		// variable

		dims.add(timeDim);
		dims.add(layerDim);
		dims.add(latDim);
		dims.add(lonDim);

		// Create variables in the output file

		Variable time = writer.addVariable(null, outTimeName, DataType.DOUBLE,
				Arrays.asList(timeDim));
		Variable layer = writer.addVariable(null, outLayerName, DataType.DOUBLE,
				Arrays.asList(layerDim));
		Variable lat = writer.addVariable(null, outLatName, DataType.DOUBLE,
				Arrays.asList(latDim));
		Variable lon = writer.addVariable(null, outLonName, DataType.DOUBLE,
				Arrays.asList(lonDim));

		// outfile.setLargeFile(true);
		var = writer.addVariable(null, outVarName, DataType.FLOAT, dims);

		// Add attribute information (cloned from source)

		cloneAttributes(uFile, inTimeName, writer, time);
		cloneAttributes(uFile, inLayerName, writer, layer);
		cloneAttributes(uFile, inLatName, writer, lat);
		cloneAttributes(uFile, inLonName, writer, lon);

		// Finalizes the structure of the output file, making the changes real.

		writer.create();

		// Write the static information for 1D axes.

		writer.write(time, uFile.findVariable(inTimeName).read());
		writer.write(layer, uFile.findVariable(inLayerName).read());
		writer.write(lat, uFile.findVariable(inLatName).read());
		writer.write(lon, uFile.findVariable(inLonName).read());

		Variable u = uFile.findVariable(inUName);
		Variable v = vFile.findVariable(inVName);

		Variable latdim = uFile.findVariable(inLatName);
		Variable londim = uFile.findVariable(inLonName);
		
		// Write the collar as NaNs.

		writeCollar(uLayerLength, u.getShape());

		// Calculate w values for each time, depth, lat and lon
		// using 3x3 (horizontal) vertical pencils

		for (int t = 0; t < u.getShape(0); t++) {
			for (int i = 1; i < u.getShape(2) - 1; i++) {
				Array lta = latdim.read(new int[] { i - 1 }, // move for
																// speed
						new int[] { 3 });
				for (int j = 1; j < u.getShape(3) - 1; j++) {

					Array lna = londim.read(new int[] { j - 1 },
							new int[] { 3 });
					
					Array uarr = u.read(new int[]{t,0,i-1,j-1}, new int[]{1,u.getShape()[1],3,3});
					Array varr = v.read(new int[]{t,0,i-1,j-1}, new int[]{1,u.getShape()[1],3,3});

					float[] w_arr = integrate(uarr, varr, lna, lta);

					Array warra = Array.factory(java.lang.Float.class,
							new int[] { 1, uLayerLength, 1, 1 }, w_arr);

					writer.write(var, new int[] { t, 0, i, j }, warra);
				}
				System.out.println("\tRow " + i + " complete.");
			}
			System.out.printf("Time %d of " + (u.getShape(0))
					+ " is complete.\n", t + 1);
		}
		return new NetcdfDataset(writer.getNetcdfFile());
	}

	
	public float[] integrate(Array u, Array v, Array lons, Array lats) {
		
		int zdim = u.getShape()[1];
		float[] w_arr = new float[zdim];
		boolean alwaysNaN = true;

		for (int k = zdim - 1; k >= 0; k--) {

			// If we're at the top layer, check if all values have
			// been NaN. If so, write NaN, otherwise, write 0.

			if (k == 0) {
				if (!alwaysNaN) {
					w_arr[k] = 0;
				} else {
					w_arr[k] = Float.NaN;
				}
				continue;
			}

			// Read the 3x3 u and v kernels

			Array ua = null;
			Array va = null;		
			
			try {
				ua = u.section(new int[]{0,k,0,0}, new int[]{1,1,3,3});
				va = v.section(new int[]{0,k,0,0}, new int[]{1,1,3,3});
			} catch (InvalidRangeException e) {
				e.printStackTrace();
			}

			float dwdz = calcdwdz(lats, lons, ua, va);

			if (Float.isNaN(dwdz)) {
				w_arr[k] = Float.NaN;
				continue;
			}

			if (alwaysNaN) {
				w_arr[k] = dwdz;
				alwaysNaN = false;
			} else {
				w_arr[k] = w_arr[k + 1] + dwdz;
			}
		}
		return w_arr;
	}

	public void go() {

		System.out.println("Writing to " + outputWFile + "...");

		try {
			uFile = NetcdfDataset.openDataset(inputUFile);
			vFile = NetcdfDataset.openDataset(inputVFile);
			generate(uFile, vFile);

		} catch (IOException e) {
			e.printStackTrace();
		} catch (InvalidRangeException e) {
			e.printStackTrace();
		}

		System.out.println("Complete.");
	}

	public void writeCollar(int layerLength, int[] shape) {

		float[] nan = new float[layerLength];
		for (int n = 0; n < layerLength; n++) {
			nan[n] = Float.NaN;
		}

		Array nana = Array.factory(java.lang.Float.class, new int[] { 1,
				layerLength, 1, 1 }, nan);
		
		try {
			for (int t = 0; t < shape[0]; t++) {
				for (int j = 0; j < shape[3]; j++) {
					writer.write(var, new int[] { t, 0, 0, j }, nana);
					writer.write(var, new int[] { t, 0, shape[2] - 1, j },
							nana);
				}
				for (int i = 1; i < shape[2] - 1; i++) {
					writer.write(var, new int[] { t, 0, i, 0 }, nana);
					writer.write(var, new int[] { t, 0, i, shape[3] - 1 },
							nana);
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		} catch (InvalidRangeException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Projects the arrays of longitudes and latitudes to a Cylindrical
	 * Equidistant Projection to bring x,y and z into a uniform coordinate
	 * system (meters).
	 * 
	 * @param lons
	 * @param lats
	 * @return
	 */

	private Array[] prj2meters(Array lons, Array lats) {
		Array[] out = new Array[2];
		Index ltidx = Index.factory(lats.getShape());
		Index lnidx = Index.factory(lons.getShape());
		Array lats_prj = Array.factory(java.lang.Double.class, lats.getShape());
		Array lons_prj = Array.factory(java.lang.Double.class, lons.getShape());
		for (int i = 0; i < lats.getSize(); i++) {
			for (int j = 0; j < lons.getSize(); j++) {
				float lon_dd = lons.getFloat(lnidx.set(j));
				float lat_dd = lats.getFloat(ltidx.set(i));
				double[] prj = MatrixUtilities.lonlat2ceqd(new double[] {
						lon_dd, lat_dd });
				lons_prj.setDouble(lnidx.set(j), prj[0]);
				lats_prj.setDouble(ltidx.set(i), prj[1]);
			}
		}
		out[0] = lons_prj;
		out[1] = lats_prj;
		return out;
	}

	/**
	 * Retrieves the name of the Latitude Variable
	 * 
	 * @return
	 */

	public String getInLatName() {
		return inLatName;
	}

	/**
	 * Retrieves the name of the Depth Variable
	 * 
	 * @return
	 */

	public String getInLayerName() {
		return inLayerName;
	}

	/**
	 * Retrieves the name of the Longitude Variable
	 * 
	 * @return
	 */

	public String getInLonName() {
		return inLonName;
	}

	/**
	 * Retrieves the u-velocity input file
	 * 
	 * @return
	 */

	public String getInputUFile() {
		return inputUFile;
	}

	/**
	 * Retrieves the v-velocity input file
	 * 
	 * @return
	 */

	public String getInputVFile() {
		return inputVFile;
	}

	/**
	 * Retrieves the name of the Time Variable
	 * 
	 * @return
	 */

	public String getInTimeName() {
		return inTimeName;
	}

	/**
	 * Retrieves the name of the u (east-west velocity) Variable
	 * 
	 * @return
	 */

	public String getInUName() {
		return inUName;
	}

	/**
	 * Retrieves the name of the v (east-west velocity) Variable
	 * 
	 * @return
	 */

	public String getInVName() {
		return inVName;
	}

	/**
	 * Retrieves the name of the output Latitude Variable
	 * 
	 * @return
	 */

	public String getOutLatName() {
		return outLatName;
	}

	/**
	 * Retrieves the name of the output Depth Variable
	 * 
	 * @return
	 */

	public String getOutLayerName() {
		return outLayerName;
	}

	/**
	 * Retrieves the name of the output Longitude Variable
	 * 
	 * @return
	 */

	public String getOutLonName() {
		return outLonName;
	}

	/**
	 * Retrieves the name of the output w file
	 * 
	 * @return
	 */

	public String getOutputWFile() {
		return outputWFile;
	}

	/**
	 * Retrieves the name of the output Time Variable
	 * 
	 * @return
	 */

	public String getOutTimeName() {
		return outTimeName;
	}

	/**
	 * Retrieves the name of the output Variable
	 * 
	 * @return
	 */

	public String getOutVarName() {
		return outVarName;
	}

	/**
	 * Sets the name of the input Latitude Variable
	 * 
	 * @param inLatName
	 */

	public void setInLatName(String inLatName) {
		this.inLatName = inLatName;
	}

	/**
	 * Sets the name of the input Depth Variable
	 * 
	 * @param inLayerName
	 */

	public void setInLayerName(String inLayerName) {
		this.inLayerName = inLayerName;
	}

	/**
	 * Sets the name of the input Longitude Variable
	 * 
	 * @param inLonName
	 */

	public void setInLonName(String inLonName) {
		this.inLonName = inLonName;
	}

	/**
	 * Sets the name of the input u (east-west) velocity file
	 * 
	 * @param inputUFile
	 */

	public void setInputUFile(String inputUFile) {
		this.inputUFile = inputUFile;
	}

	/**
	 * Sets the name of the input v (north-south) velocity file
	 * 
	 * @param inputVFile
	 */

	public void setInputVFile(String inputVFile) {
		this.inputVFile = inputVFile;
	}

	/**
	 * Sets the name of the input Time Variable
	 * 
	 * @param inTimeName
	 */

	public void setInTimeName(String inTimeName) {
		this.inTimeName = inTimeName;
	}

	/**
	 * Sets the name of the input u (east-west) velocity Variable
	 * 
	 * @param inUName
	 */

	public void setInUName(String inUName) {
		this.inUName = inUName;
	}

	/**
	 * Sets the name of the input v (north-south) velocity Variable
	 * 
	 * @param inVName
	 */

	public void setInVName(String inVName) {
		this.inVName = inVName;
	}

	/**
	 * Sets the name of the output Latitude Variable
	 * 
	 * @param outLatName
	 */

	public void setOutLatName(String outLatName) {
		this.outLatName = outLatName;
	}

	/**
	 * Sets the name of the output Depth Variable
	 * 
	 * @param outLayerName
	 */

	public void setOutLayerName(String outLayerName) {
		this.outLayerName = outLayerName;
	}

	/**
	 * Sets the name of the output Longitude Variable
	 * 
	 * @param outLonName
	 */

	public void setOutLonName(String outLonName) {
		this.outLonName = outLonName;
	}

	/**
	 * Sets the name of the output w velocity file
	 * 
	 * @param outputWFile
	 */

	public void setOutputWFile(String outputWFile) {
		this.outputWFile = outputWFile;
	}

	/**
	 * Sets the name of the output Time Variable
	 * 
	 * @param outTimeName
	 */

	public void setOutTimeName(String outTimeName) {
		this.outTimeName = outTimeName;
	}

	/**
	 * Sets the name of the output Variable
	 * 
	 * @param outVarName
	 */

	public void setOutVarName(String outVarName) {
		this.outVarName = outVarName;
	}

	public void setReproject(boolean reproject) {
		this.reproject = reproject;
	}
}