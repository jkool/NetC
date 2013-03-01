package conv;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import ucar.ma2.Array;
import ucar.ma2.DataType;
import ucar.ma2.InvalidRangeException;
import ucar.nc2.Attribute;
import ucar.nc2.Dimension;
import ucar.nc2.NetcdfFile;
import ucar.nc2.NetcdfFileWriteable;
import ucar.nc2.Variable;
import ucar.nc2.dataset.NetcdfDataset;

public class NCFBlock {

	private double maxFileSize = 4E9;
	private String fileDir = "V:/HYCOM";
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
	private String outPath = "C:/Temp/testfile.nc";
	private NetcdfFile ncf_in;
	private Array timeArr;
	private Array depthArr;
	private Array latArr;
	private Array lonArr;
	private int[] shape;
	private boolean negativeDepth = true;

	public static void main(String[] args) {
		NCFBlock ncb = new NCFBlock();
		try {
			ncb.readTemplate("V:/HYCOM/AUS_w_2009_1.nc");
			NetcdfFileWriteable nc = ncb.createOutput("Y:/NERP_metadata/Samples/w_velocity.nc");
			ncb.writeOutput(nc);
			
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.out.println("Complete.");
	}

	public void readTemplate(String inFileName) {
		try {
			ncf_in = NetcdfFile.open(inFileName);
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
				System.out.println("WARNING: Depth variable " + inDepthName
						+ " not found, and the number of dimensions is greater than 3."
						+ "  File variables are:\n"
						+ Arrays.toString(var.toArray()));
			}

			int latidx = lat.findDimensionIndex(inVerticalDim);
			int lonidx = lat.findDimensionIndex(inHorizontalDim);
			int latlen = lat.getDimension(latidx).getLength();
			int lonlen = lon.getDimension(lonidx).getLength();

			timeArr = time.read();
			
			if(depth!=null){
				depthArr = depth.read();
				if(negativeDepth){
					for(int i = 0; i < depthArr.getShape()[0]; i++){
						if(i!=0){
						depthArr.setDouble(i, -depthArr.getDouble(i));
						}else{depthArr.setDouble(i, 0);}
					}
				}
			}

			int[] latshape = ones(lat.getRank());
			latshape[latidx] = latlen;
			latArr = (lat.read(new int[lat.getRank()], latshape)).reduce();

			int[] lonshape = ones(lon.getRank());
			lonshape[lonidx] = lonlen;
			lonArr = (lon.read(new int[lon.getRank()], lonshape)).reduce();

		} catch (IOException e) {
			e.printStackTrace();
		} catch (InvalidRangeException e) {
			e.printStackTrace();
		}
	}

	private NetcdfFileWriteable createOutput(String outputFile) {
		NetcdfFileWriteable ncf_out = null;
		try {
			ncf_out = NetcdfFileWriteable.createNew(outputFile,false);
			
			// Add Dimensions
			Dimension timeDim = new Dimension(outTimeName, timeArr.getShape()[0]);
			Dimension depthDim = null;
			if(depthArr!=null){
				depthDim = new Dimension(outDepthName, depthArr.getShape()[0]);
			}
			
			Dimension latDim = new Dimension(outLatName, latArr.getShape()[0]);
			Dimension lonDim = new Dimension(outLonName, lonArr.getShape()[0]);
			
			ncf_out.addDimension(null, timeDim);
			if(depthDim!=null){
				ncf_out.addDimension(null,  depthDim);
			}
			ncf_out.addDimension(null,  latDim);
			ncf_out.addDimension(null,  lonDim);
			
			// Add Variables
			ncf_out.addVariable(outTimeName, DataType.DOUBLE, new Dimension[]{timeDim});
			
			if(depthDim!=null){
				ncf_out.addVariable(outDepthName, DataType.DOUBLE, new Dimension[]{depthDim});
			}
			
			ncf_out.addVariable(outLatName, DataType.DOUBLE, new Dimension[]{latDim});
			ncf_out.addVariable(outLonName, DataType.DOUBLE, new Dimension[]{lonDim});
			
			Dimension[] dims = null;
			if(depthDim==null){
				dims = new Dimension[]{timeDim,latDim,lonDim};
			}
			else{
				dims = new Dimension[]{timeDim,depthDim,latDim,lonDim};
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
		}
		return ncf_out;
	}
	
	private void writeOutput(NetcdfFileWriteable ncf_out){
		try {

			// Write the static information for 1D axes.

			ncf_out.write(outTimeName, timeArr);
			ncf_out.write(outDepthName, depthArr);
			ncf_out.write(outLatName, latArr);
			ncf_out.write(outLonName, lonArr);
			
			int timesteps = ncf_out.findDimension(outTimeName).getLength();
			Variable invar = ncf_in.findVariable(inVarName);
			Variable outvar = ncf_out.findVariable(outVarName);
			
			int[] shape = outvar.getShape();
			shape[0] = 1;
			
			for(int i = 0; i < timesteps; i++){
				int[] origin = new int[]{i,0,0,0};
				ncf_out.write(outVarName, origin, invar.read(origin,shape));
			}
			
			ncf_out.flush();
			ncf_out.close();
			
			
			
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (InvalidRangeException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
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

	private int[] ones(int length) {
		int[] out = new int[length];
		for (int i = 0; i < length; i++) {
			out[i] = 1;
		}
		return out;
	}
}
