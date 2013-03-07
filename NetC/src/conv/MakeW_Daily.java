package conv;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import ucar.ma2.*;
import ucar.nc2.*;
import ucar.nc2.dataset.NetcdfDataset;
import utilities.MatrixUtilities;

/**
 * Converts daily u and v netCDF velocity files into
 * w velocity files.  The Daily routine is different in that
 * Latitude and Longitude are 2 dimensional, and the dimensions are
 * based on a X and Y index system.
 */

public class MakeW_Daily {

	public static void main(String[] args) {
		MakeW_Daily mw = new MakeW_Daily();
		mw.go();
	}

	private String inputUFile = "V:/HYCOM/AUS_u_2009_4.nc";;
	private String inputVFile = "V:/HYCOM/AUS_v_2009_4.nc";;
	private String outputWFile = "V:/HYCOM/AUS_w_2009_4.nc";
	private String inputBathyFile = "E:/HPC/Modeling/AUS/Input/NetCDF/GLB_depth_2d.nc";
	private ArrayList<Dimension> dims;
	private NetcdfFileWriteable outfile;
	private NetcdfDataset uFile, vFile;
	private NetCDF_FloatGrid_3D bathy;
	private String inLatName = "Latitude";
	private String inLonName = "Longitude";
	private String inLayerName = "Depth";
	private String inTimeName = "MT";
	private String inXName = "X";
	private String inYName = "Y";
	private String inUName = "u";
	private String inVName = "v";
	private String outLatName = inLatName;
	private String outLonName = inLonName;
	private String outLayerName = inLayerName;
	private String outTimeName = inTimeName;
	private String outXName = inXName;
	private String outYName = inYName;
	private String outVarName = "w";
	private float bathy_cutoff = 1E30f;
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

	private float calcdwdz(Array lats, Array lons, Array us, Array vs) {
		
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
		
		double dx0 = prj[0].getDouble(0);
		double dx1 = prj[0].getDouble((int) prj[0].getSize() - 1);
		double dx = dx1 - dx0;
		double dy0 = prj[1].getDouble(1);
		double dy1 = prj[1].getDouble((int) prj[1].getSize() - 1);
		double dy = dy1 - dy0;
		float dudx = du / (float) dx;
		float dvdy = dv / (float) dy;
		return dudx + dvdy;
	}
	
	/**
	 * Calculates the difference in depth
	 * 
	 * @param k - the index of the depth layer
	 * @param mpz - midpoint depth values
	 * @param maxZ - maximum depth value
	 * @return
	 */

	private float calcdz(int k, float[] mpz, float maxZ) {
		if (k == mpz.length) {
			return maxZ - mpz[k - 1];
		}

		if (maxZ < mpz[k]) {
			return maxZ - mpz[k - 1];
		} else {
			return mpz[k] - mpz[k - 1];
		}
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
			NetcdfFileWriteable outfile, String outVarName) {

		// Find the variable

		Variable vi = infile.findVariable(inVarName);

		List<Attribute> l = vi.getAttributes();

		for (Attribute a : l) {
			outfile.addVariableAttribute(outVarName, a);
		}
	}

	/**
	 * Calculates change in the x direction
	 * 
	 * @param arr
	 *            - a 3x3 array of numbers
	 * @return
	 */

	private float dx(Array arr) {
		Index idx = Index.factory(arr.getShape());
		float e = arr.getFloat(idx.set(1, 1));
		if (Float.isNaN(e)) {
			return Float.NaN;
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
	 * Calculates change in the y direction
	 * 
	 * @param arr
	 *            - a 3x3 array of numbers
	 * @return
	 */

	private float dy(Array arr) {
		Index idx = Index.factory(arr.getShape());
		float e = arr.getFloat(idx.set(1, 1));
		if (Float.isNaN(e)) {
			return Float.NaN;
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

	private void generate(NetcdfDataset uFile, NetcdfDataset vFile)
			throws InvalidRangeException, IOException {

		// Set up empty arrays for ranges and dimensions.

		dims = new ArrayList<Dimension>();

		// Create the output file

		outfile = NetcdfFileWriteable.createNew(outputWFile, false);

		// Construct the data set dimensions - Time, Depth, Latitude and
		// Longitude (in order)

		int uTimeLength = uFile.findDimension(inTimeName).getLength();
		int uLayerLength = uFile.findDimension(inLayerName).getLength();
		int uXLength = uFile.findDimension(inXName).getLength();
		int uYLength = uFile.findDimension(inYName).getLength();

		int vTimeLength = vFile.findDimension(inTimeName).getLength();
		int vLayerLength = vFile.findDimension(inLayerName).getLength();
		int vXLength = vFile.findDimension(inXName).getLength();
		int vYLength = vFile.findDimension(inYName).getLength();

		Dimension timeDim = outfile.addDimension(outTimeName,
				Math.min(uTimeLength, vTimeLength));
		Dimension layerDim = outfile.addDimension(outLayerName,
				Math.min(uLayerLength, vLayerLength));
		Dimension yDim = outfile.addDimension(outYName,
				Math.min(uYLength, vYLength));
		Dimension xDim = outfile.addDimension(outXName,
				Math.min(uXLength, vXLength));

		// Add to a list - this becomes the coordinate system for the output
		// variable. The order is *very* important!

		dims.add(timeDim);
		dims.add(layerDim);
		dims.add(yDim);
		dims.add(xDim);

		// Create variables in the output file

		outfile.addVariable(outTimeName, DataType.DOUBLE,
				new Dimension[] { timeDim });
		outfile.addVariable(outLayerName, DataType.DOUBLE,
				new Dimension[] { layerDim });
		outfile.addVariable(outLatName, DataType.DOUBLE, new Dimension[] {
				yDim, xDim });
		outfile.addVariable(outLonName, DataType.DOUBLE, new Dimension[] {
				yDim, xDim });

		// outfile.setLargeFile(true);
		outfile.addVariable(outVarName, DataType.FLOAT, dims);

		// Add attribute information (cloned from source)

		cloneAttributes(uFile, inTimeName, outfile, outTimeName);
		cloneAttributes(uFile, inLayerName, outfile, outLayerName);
		cloneAttributes(uFile, inLatName, outfile, outLatName);
		cloneAttributes(uFile, inLonName, outfile, outLonName);

		// Finalizes the structure of the output file, making the changes real.

		outfile.create();

		// Write the static information for 1D axes.

		outfile.write(outTimeName, uFile.findVariable(inTimeName).read());
		outfile.write(outLayerName, uFile.findVariable(inLayerName).read());
		outfile.write(outLatName, uFile.findVariable(inLatName).read());
		outfile.write(outLonName, uFile.findVariable(inLonName).read());

		Variable u = uFile.findVariable(inUName);
		Variable v = vFile.findVariable(inVName);

		Variable latdim = uFile.findVariable(inLatName);
		Variable londim = uFile.findVariable(inLonName);
		Variable zdim = uFile.findVariable(inLayerName);

		// Read the depth array

		Array depth_array = zdim.read(new int[] { 0 },
				new int[] { uLayerLength });
		
		double[] depths;

		if (depth_array.getClass().getName().contains("Float")) {

			float[] tmp = (float[]) depth_array.copyTo1DJavaArray();
			depths = new double[tmp.length];

			for (int i = 0; i < tmp.length; i++) {
				depths[i] = tmp[i];
			}
		}

		else {
			depths = (double[]) depth_array.copyTo1DJavaArray();
		}

		// Determine the size and midpoint of the depth bins

		double[] binz = new double[depths.length - 1];
		float[] mpz = new float[depths.length - 1];

		for (int i = 0; i < depths.length - 1; i++) {
			binz[i] = depths[i + 1] - depths[i];
			mpz[i] = (float) (depths[i + 1] + depths[i]) / 2;
		}

		float[] nan = new float[uLayerLength];
		for (int n = 0; n < uLayerLength; n++) {
			nan[n] = Float.NaN;
		}

		// Write the collar as NaNs.

		Array nana = Array.factory(java.lang.Float.class, new int[] { 1,
				uLayerLength, 1, 1 }, nan);
		for (int t = 0; t < u.getShape(0); t++) {
			for (int j = 0; j < u.getShape(3); j++) {
				outfile.write(outVarName, new int[] { t, 0, 0, j }, nana);
				outfile.write(outVarName, new int[] { t, 0, u.getShape(2) - 1,
						j }, nana);
			}
			for (int i = 1; i < u.getShape(2) - 1; i++) {
				outfile.write(outVarName, new int[] { t, 0, i, 0 }, nana);
				outfile.write(outVarName, new int[] { t, 0, i,
						u.getShape(3) - 1 }, nana);
			}
		}

		// Calculate w values for each time, depth, lat and lon
		// using 3x3 (horizontal) vertical pencils

		for (int t = 0; t < u.getShape(0); t++) {
			for (int i = 1; i < u.getShape(2) - 1; i++) {
				Array lta = latdim.read(new int[] { i - 1, 0 }, // move for
																// speed
						new int[] { 3, 1 }).reduce();
				double currentLat = lta.getDouble(1);
				for (int j = 1; j < u.getShape(3) - 1; j++) {

					float[] w_arr = new float[uLayerLength];

					if (i == 0 || j == 0 || i == u.getShape(2) - 1
							|| j == u.getShape(3) - 1) {

						for (int n = 0; n < w_arr.length; n++) {
							w_arr[n] = Float.NaN;
						}
						continue;
					}

					Array lna = londim.read(new int[] { 0, j - 1 },
							new int[] { 1, 3 }).reduce();

					double currentLon = lna.getDouble(1);
					float maxZ = bathy.value_at_lonlat((float) currentLon,
							(float) currentLat);

					boolean alwaysNaN = true;

					for (int k = u.getShape(1) - 1; k >= 0; k--) {

						// If the bathymetry says we're on land, then all values
						// will be NaN. This type of check is required because
						// the value does not equal java.Float.NaN, but
						// 1.2676506E30.

						if (maxZ > bathy_cutoff) {
							w_arr[k] = Float.NaN;
							continue;
						}

						if (k == 0) {
							if (!alwaysNaN) {
								w_arr[k] = 0;
							} else {
								w_arr[k] = Float.NaN;
							}
							continue;
						}

						if (maxZ < mpz[k - 1]) {
							w_arr[k] = Float.NaN;
							continue;
						}

						Array ua = u.read(new int[] { t, k, i - 1, j - 1 },
								new int[] { 1, 1, 3, 3 }).reduce();
						Array va = v.read(new int[] { t, k, i - 1, j - 1 },
								new int[] { 1, 1, 3, 3 }).reduce();

						float dz = calcdz(k, mpz, maxZ);

						if (Float.isNaN(dz)) {
							w_arr[k] = dz;
							continue;
						}

						float dwdz = calcdwdz(lta, lna, ua, va);
						if (Float.isNaN(dwdz)) {
							continue;
						}

						if (alwaysNaN) {
							w_arr[k] = dwdz;
							alwaysNaN = false;
						} else {
							w_arr[k] = w_arr[k + 1] + dwdz;
						}
					}
					Array warra = Array.factory(java.lang.Float.class,
							new int[] { 1, uLayerLength, 1, 1 }, w_arr);

					outfile.write(outVarName, new int[] { t, 0, i, j }, warra);
				}
				System.out.println("\tRow " + i + " complete.");
			}
			System.out.printf("Time %d of " + (u.getShape(0))
					+ " is complete.\n", t + 1);
		}
	}

	/**
	 * Main-type execution
	 */

	public void go() {

		System.out.println("Writing to " + outputWFile + "...");

		try {
			uFile = NetcdfDataset.openDataset(inputUFile);
			vFile = NetcdfDataset.openDataset(inputVFile);
			bathy = new NetCDF_FloatGrid_3D(inputBathyFile, inLatName,
					inLonName);
			generate(uFile, vFile);

		} catch (IOException e) {
			e.printStackTrace();
		} catch (InvalidRangeException e) {
			e.printStackTrace();
		}

		System.out.println("Complete.");
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
		assert Arrays.equals(lats.getShape(), lons.getShape());
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
	 * Returns the value used as a cutoff for determining valid bathymetry
	 * values.
	 * 
	 * @return
	 */

	public float getBathy_cutoff() {
		return bathy_cutoff;
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
	 * Retrieves the name of the bathymetry file used as a boundary.
	 * 
	 * @return
	 */

	public String getInputBathyFile() {
		return inputBathyFile;
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
	 * Retrieves the name of the horizontal Variable
	 * 
	 * @return
	 */

	public String getInXName() {
		return inXName;
	}

	/**
	 * Retrieves the name of the vertical (2D) Variable
	 * 
	 * @return
	 */

	public String getInYName() {
		return inYName;
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
	 * Retrieves the name of the output horizontal Variable
	 * 
	 * @return
	 */

	public String getOutXName() {
		return outXName;
	}

	/**
	 * Retrieves the name of the output vertical (2D) Variable
	 * 
	 * @return
	 */

	public String getOutYName() {
		return outYName;
	}
	
	/**
	 * Sets the bathymetry cutoff value
	 * 
	 * @param bathy_cutoff
	 *            - value beyond which bathymetry values are considered to be
	 *            non-values
	 */

	public void setBathy_cutoff(float bathy_cutoff) {
		this.bathy_cutoff = bathy_cutoff;
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
	 * Sets the name of the input bathymetry file
	 * 
	 * @param inputBathyFile
	 */

	public void setInputBathyFile(String inputBathyFile) {
		this.inputBathyFile = inputBathyFile;
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
	 * Sets the name of the input horizontal Variable
	 * 
	 * @param inXName
	 */

	public void setInXName(String inXName) {
		this.inXName = inXName;
	}

	/**
	 * Sets the name of the input vertical (2D) Variable
	 * 
	 * @param inYName
	 */

	public void setInYName(String inYName) {
		this.inYName = inYName;
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

	/**
	 * Sets the name of the output horizontal Variable
	 * 
	 * @param outXName
	 */

	public void setOutXName(String outXName) {
		this.outXName = outXName;
	}

	/**
	 * Sets the name of the output vertical(2D) variable
	 * 
	 * @param outYName
	 */

	public void setOutYName(String outYName) {
		this.outYName = outYName;
	}
	
	/**
	 * Sets whether the data should be re-projected when
	 * calculating the velocity values
	 * 
	 * @param reproject
	 */
	
	public void setReproject(boolean reproject){
		this.reproject = reproject;
	}
}
