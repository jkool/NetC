package conv;

import java.io.IOException;
import java.util.Arrays;

import ucar.ma2.Array;
import ucar.ma2.InvalidRangeException;
import ucar.nc2.NetcdfFile;
import ucar.nc2.Variable;

public class NetCDF_FloatGrid_3D {

	private NetcdfFile file;
	private Variable var;
	private Variable latvar;
	private Variable lonvar;
	private String variableName = "bathymetry";
	private String latName = "Latitude";
	private String lonName = "Longitude";
	private float[] lats, lons;
	private Array array, latarr, lonarr;

	public NetCDF_FloatGrid_3D(String fileName) {

		try {
			file = NetcdfFile.open(fileName);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		var = file.findVariable(variableName);
		latvar = file.findVariable(latName);
		lonvar = file.findVariable(lonName);
		
		int[] shp = latvar.getShape();
		
		try {
			if (latvar.getRank() == 2) {
				latarr = latvar.read(new int[] { 0, 0 },
						new int[] { shp[0], 1 });
			} else {
				latarr = latvar.read(new int[] { 0 }, new int[] { shp[0] });
			}
			lats = toJArray(latarr);

			if (lonvar.getRank() == 2) {
				lonarr = lonvar.read(new int[] { 0, 0 }, new int[] { 1,
						shp[1] });
			} else {
				lonarr = lonvar.read(new int[] { 0 },
						new int[] { lonvar.getShape()[0] });
			}
			lons = toJArray(lonarr);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (InvalidRangeException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public NetCDF_FloatGrid_3D(String fileName, String latName, String lonName) {

		try {
			file = NetcdfFile.open(fileName);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		this.latName = latName;
		this.lonName = lonName;
		var = file.findVariable(variableName);
		latvar = file.findVariable(latName);
		lonvar = file.findVariable(lonName);
		int[] shp = latvar.getShape();
		
		try {
			if (latvar.getRank() == 2) {
				latarr = latvar.read(new int[] { 0, 0 },
						new int[] { shp[0], 1 });
			} else {
				latarr = latvar.read(new int[] { 0 }, new int[] { shp[0] });
			}
			lats = toJArray(latarr);

			if (lonvar.getRank() == 2) {
				lonarr = lonvar.read(new int[] { 0, 0 }, new int[] { 1,
						shp[1] });
			} else {
				lonarr = lonvar.read(new int[] { 0 },
						new int[] { lonvar.getShape()[0] });
			}
			lons = toJArray(lonarr);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (InvalidRangeException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	// With single dimension lat and lon, this only needs to be done once.!!!!!  ALso check for land zeroes and convert to NaN.
	
	public int[] lonlat2ij(float lon, float lat) {

		int[] result = new int[] { -1, -1 };

		result[0] = binary_search(lats, lat);	
		result[1] = binary_search(lons, lon < lons[0] ? lon + 360 : lon);

		return result;
	}

	public float value_at_ij(int i, int j) {
		try {
			if (var.getRank() == 2) {
				array = var.read(new int[] { i, j }, new int[] { 1, 1 });
			} else {
				array = var.read(new int[] { 0, i, j }, new int[] { 1, 1, 1 });
			}
			return array.getFloat(0);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (InvalidRangeException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return Float.NaN;
	}

	public float value_at_lonlat(float lon, float lat) {
		int[] idx = lonlat2ij(lon, lat);
		return value_at_ij(idx[0], idx[1]);
	}

	private float[] toJArray(Array arr) {
		float[] ja;
		if (arr.getElementType() == Double.TYPE) {

			double[] da = (double[]) arr.copyTo1DJavaArray();
			ja = new float[da.length];
			for (int i = 0; i < ja.length; i++) {
				ja[i] = (float) da[i];
			}
			return ja;
		} else {

			return (float[]) arr.copyTo1DJavaArray();
		}
	}

	private int binary_search(float[] array, float val) {
		int idx = Arrays.binarySearch(array, val);

		if (idx < 0) {

			// Error check

			if (idx == -1) {
				throw new IllegalArgumentException(var.getShortName()
						+ " value " + val + " does not fall in the range "
						+ array[0] + " : " + array[array.length - 1] + ".");
			}

			return -(idx + 2);
		}

		return idx;
	}

	public String getVariableName() {
		return variableName;
	}

	public void setVariableName(String variableName) {
		this.variableName = variableName;
	}
}
