package conv;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import ucar.ma2.Array;
import ucar.ma2.DataType;
import ucar.ma2.InvalidRangeException;
import ucar.nc2.Attribute;
import ucar.nc2.Dimension;
import ucar.nc2.NetcdfFile;
import ucar.nc2.NetcdfFileWriteable;
import ucar.nc2.Variable;

public class NetCDFDir2Block {

	private double maxFileSize = 4E9;
	private String fileDir = "V:/Data/HYCOM";
	private String inVarName = "w";
	private String inHorizontalDim = "X";
	private String inVerticalDim = "Y";
	private String inTimeName = "MT";
	private String inDepthName = "Depth";
	private String inLatName = "Latitude";
	private String inLonName = "Longitude";
	private String outTimeName = "Time";
	private String outDepthName = "Depth";
	private String outLatName = inLatName;
	private String outLonName = inLonName;
	private String outVarName = inVarName;
	private String outputDir = "V:/Data/HYCOM/Blocks";
	private DateFormat df = new SimpleDateFormat("yyyy_MM_dd");
	private FilenameFilter filter = new FilenamePatternFilter(".*_["
			+ inVarName + "]_.*\\.nc");
	private long blocksize = TimeConvert.daysToMillis("30");
	private Array timeArr;
	private Array depthArr;
	private Array latArr;
	private Array lonArr;
	private boolean negativeDepth = true;
	private NetcdfFile ncf_in;

	public static void main(String[] args) {
		NetCDFDir2Block ncb = new NetCDFDir2Block();
		ncb.convert(ncb.fileDir, ncb.filter, ncb.blocksize);
		System.out.println("Complete.");
	}

	public void convert(String fileDir, FilenameFilter filter, long blocksize) {
		this.fileDir = fileDir;
		this.filter = filter;
		this.blocksize = blocksize;
		convert();
	}

	public void convert() {
		NetCDFDir ncdir = new NetCDFDir(fileDir, filter);
		ArrayList<Long> times = new ArrayList<Long>(ncdir.getTimes());
		List<Long> blocks = calcBlock(times, blocksize);

		TreeMap<Long, NetcdfFile> selection;

		// iterate through the blocks
		for (int i = 1; i < blocks.size(); i++) {
			System.out.println("Working on block " + i + "...");
			long start_time = blocks.get(i - 1);
			long end_time = blocks.get(i);
			selection = selectBetween(start_time, end_time, ncdir.getFiles());
			String outputPath = outputDir + "/"
					+ df.format(new Date(start_time)) + "_to_"
					+ df.format(new Date(end_time)) + "_tmp_" + inVarName + ".nc";
			NetcdfFileWriteable ncf_write = makeTemplate(selection, outputPath);
			copyContent(selection, ncf_write);
			close(ncf_write);
			System.out.println(outputPath + " has been written.");
		}

		// get the remainder
		long time = blocks.get(blocks.size() - 1);
		selection = selectAfter(time, ncdir.getFiles());
		String outputPath = outputDir + "/" + df.format(new Date(time))
				+ "_onward_" + inVarName + ".nc";
		NetcdfFileWriteable ncf_write = makeTemplate(selection, outputPath);
		copyContent(selection, ncf_write);
	}

	public NetcdfFileWriteable makeTemplate(TreeMap<Long, NetcdfFile> files,
			String outputFile) {

		NetcdfFileWriteable ncf_out = null;
		if (!(new File(outputDir).exists())) {
			throw new IllegalArgumentException("Directory " + outputDir
					+ " does not exist.");
		}

		try {
			ncf_in = files.get((new ArrayList<Long>(files.keySet())).get(0));

			Variable time = ncf_in.findVariable(inTimeName);
			if (time == null) {
				List<Variable> var = ncf_in.getVariables();
				System.out.println("WARNING: Variable " + inTimeName
						+ " was not found.  File variables are:\n"
						+ Arrays.toString(var.toArray()));
			}

			Variable depth = ncf_in.findVariable(inDepthName);

			if (depth == null && ncf_in.getDimensions().size() > 3) {
				List<Variable> var = ncf_in.getVariables();
				System.out
						.println("WARNING: Depth variable "
								+ inDepthName
								+ " not found, and the number of dimensions is greater than 3."
								+ "  File variables are:\n"
								+ Arrays.toString(var.toArray()));
			}

			Variable lat = ncf_in.findVariable(inLatName);
			if (lat == null) {
				List<Variable> var = ncf_in.getVariables();
				System.out.println("WARNING: Variable " + inLatName
						+ " was not found.  File variables are:\n"
						+ Arrays.toString(var.toArray()));
			}

			Variable lon = ncf_in.findVariable(inLonName);
			if (lon == null) {
				List<Variable> var = ncf_in.getVariables();
				System.out.println("WARNING: Variable " + inLonName
						+ " was not found.  File variables are:\n"
						+ Arrays.toString(var.toArray()));
			}

			int latidx = lat.findDimensionIndex(inVerticalDim);
			int lonidx = lat.findDimensionIndex(inHorizontalDim);
			int latlen = lat.getDimension(latidx).getLength();
			int lonlen = lon.getDimension(lonidx).getLength();

			Long[] la = files.keySet().toArray(new Long[files.size()]);
			timeArr = Array
					.factory(DataType.DOUBLE, new int[] { files.size() });

			for (int i = 0; i < la.length; i++) {
				timeArr.setDouble(i, TimeConvert.millisToHYCOM(la[i]));
			}

			if (depth != null) {
				depthArr = depth.read();
				if (negativeDepth) {
					for (int i = 0; i < depthArr.getShape()[0]; i++) {
						if (i != 0) {
							depthArr.setDouble(i, -depthArr.getDouble(i));
						} else {
							depthArr.setDouble(i, 0);
						}
					}
				}
			}

			int[] latshape = ones(lat.getRank());
			latshape[latidx] = latlen;
			latArr = (lat.read(new int[lat.getRank()], latshape)).reduce();

			int[] lonshape = ones(lon.getRank());
			lonshape[lonidx] = lonlen;
			lonArr = (lon.read(new int[lon.getRank()], lonshape)).reduce();

			ncf_out = NetcdfFileWriteable.createNew(outputFile, false);

			// Add Dimensions
			Dimension timeDim = new Dimension(outTimeName,
					timeArr.getShape()[0]);
			Dimension depthDim = null;

			if (depthArr != null) {
				depthDim = new Dimension(outDepthName, depthArr.getShape()[0]);
			}

			Dimension latDim = new Dimension(outLatName, latArr.getShape()[0]);
			Dimension lonDim = new Dimension(outLonName, lonArr.getShape()[0]);

			ncf_out.addDimension(null, timeDim);
			if (depthDim != null) {
				ncf_out.addDimension(null, depthDim);
			}
			ncf_out.addDimension(null, latDim);
			ncf_out.addDimension(null, lonDim);

			// Add Variables
			ncf_out.addVariable(outTimeName, DataType.DOUBLE,
					new Dimension[] { timeDim });

			if (depthDim != null) {
				ncf_out.addVariable(outDepthName, DataType.DOUBLE,
						new Dimension[] { depthDim });
			}

			ncf_out.addVariable(outLatName, DataType.DOUBLE,
					new Dimension[] { latDim });
			ncf_out.addVariable(outLonName, DataType.DOUBLE,
					new Dimension[] { lonDim });

			Dimension[] dims = null;

			if (depthDim == null) {
				dims = new Dimension[] { timeDim, latDim, lonDim };
			} else {
				dims = new Dimension[] { timeDim, depthDim, latDim, lonDim };
			}

			ncf_out.addVariable(outVarName, DataType.FLOAT, dims);

			// Add attribute information (cloned from source)

			cloneAttributes(ncf_in, inTimeName, ncf_out, outTimeName);
			cloneAttributes(ncf_in, inDepthName, ncf_out, outDepthName);
			cloneAttributes(ncf_in, inLatName, ncf_out, outLatName);
			cloneAttributes(ncf_in, inLonName, ncf_out, outLonName);
			cloneAttributes(ncf_in, inVarName, ncf_out, outVarName);

			ncf_out.create();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (InvalidRangeException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		return ncf_out;
	}

	private void copyContent(TreeMap<Long, NetcdfFile> map,
			NetcdfFileWriteable ncf_out) {
		try {

			// Write the static information for 1D axes.

			ncf_out.write(outTimeName, timeArr);
			ncf_out.write(outDepthName, depthArr);
			ncf_out.write(outLatName, latArr);
			ncf_out.write(outLonName, lonArr);

			Iterator<Long> it = map.keySet().iterator();
			int ct = 0;

			while (it.hasNext()) {
				long time = it.next();
				System.out.print("\t" + ct + "\tWorking on " + new Date(time));
				NetcdfFile copyfile = map.get(time);
				int timesteps = copyfile.findDimension(inTimeName).getLength();
				Variable invar = copyfile.findVariable(inVarName);

				int[] shape = invar.getShape();
				shape[0] = 1;

				for (int i = 0; i < timesteps; i++) {
					int[] origin_in = new int[] { i, 0, 0, 0 };
					int[] origin_out = new int[] { ct, 0, 0, 0 };
					ncf_out.write(outVarName, origin_out,
							invar.read(origin_in, shape));
					ct++;
				}
				System.out.println("\tDone.");

			}

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (InvalidRangeException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private void close(NetcdfFileWriteable ncf_out) {
		if (ncf_out != null) {
			try {
				ncf_out.flush();
				ncf_out.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}

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

	private List<Long> calcBlock(List<Long> times, long duration) {
		if (times.size() == 0) {
			return null;
		}
		ArrayList<Long> blockstart = new ArrayList<Long>();
		int n = times.size();
		long bstart = times.get(0);
		blockstart.add(bstart);
		for (int i = 0; i < n; i++) {
			if (!(times.get(i) - bstart < duration)) {
				bstart = times.get(i);
				blockstart.add(times.get(i));
			}
		}
		return blockstart;
	}

	private TreeMap<Long, NetcdfFile> selectBetween(long from, long until,
			Map<Long, NetcdfFile> map) {
		TreeMap<Long, NetcdfFile> outmap = new TreeMap<Long, NetcdfFile>();
		Iterator<Long> it = map.keySet().iterator();
		while (it.hasNext()) {
			long time = it.next();
			if (time >= from && time < until) {
				outmap.put(time, map.get(time));
			}
		}
		return outmap;
	}

	private TreeMap<Long, NetcdfFile> selectAfter(long from,
			Map<Long, NetcdfFile> map) {
		TreeMap<Long, NetcdfFile> outmap = new TreeMap<Long, NetcdfFile>();
		Iterator<Long> it = map.keySet().iterator();
		while (it.hasNext()) {
			long time = it.next();
			if (time >= from) {
				outmap.put(time, map.get(time));
			}
		}
		return outmap;
	}

	private int[] ones(int length) {
		int[] out = new int[length];
		for (int i = 0; i < length; i++) {
			out[i] = 1;
		}
		return out;
	}
}
