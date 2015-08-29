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
import ucar.nc2.NetcdfFile;
import ucar.nc2.NetcdfFileWriter;
import ucar.nc2.Variable;

public class SquareCell4D {
	
	private NetcdfFile ncf;
	private final float NODATA = 1E34f;
	
	private String inLatName = "Latitude";
	private String inLonName = "Longitude";
	private String inLayerName = "Depth";
	private String inTimeName = "Time";
	private String inVarName = "salinity";
	private String outLatName = inLatName;
	private String outLonName = inLonName;
	private String outLayerName = inLayerName;
	private String outTimeName = inTimeName;
	private String outVarName = inVarName;
	private static String inFile = "Y:/NERP_metadata/Samples/salinity.nc";
	private String outFile = "Y:/NERP_metadata/Samples/salinity_sq.nc";
	private NetcdfFileWriter writer;
	private Variable outTime,outLayer,outLat,outLon,outVar;
	private float lb_valid = -300f;
	private float ub_valid = 300f;
	private float minlat,maxlat,minlon,maxlon;
	
	public static void main(String[] args){
		
		SquareCell4D sq = new SquareCell4D();
		try {
			sq.square(inFile);
		} catch (InvalidRangeException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		System.out.println("Complete.");
		System.exit(0);
	}
	
	private void square(String filename) throws InvalidRangeException, IOException{
		
		int nlat,nlon;
		
		ncf = NetcdfFile.open(filename);
		
		Variable inLat = ncf.findVariable(inLatName);
		Variable inLon = ncf.findVariable(inLonName);
		Variable inDepth = ncf.findVariable(inLayerName);
		Variable inTime = ncf.findVariable(inTimeName);
		Variable inVar = ncf.findVariable(inVarName);
	
		int dlen = inDepth.getShape()[0];
		int tlen = inTime.getShape()[0];
		
		Array latarr = inLat.read();
		Array lonarr = inLon.read();
		
		double[] platarr = (double[]) latarr.copyTo1DJavaArray();
		double[] plonarr = (double[]) lonarr.copyTo1DJavaArray();
		
		Index latindex = latarr.getIndex();
		Index lonindex = lonarr.getIndex();
		
		nlat = latarr.getShape()[0];
		nlon = lonarr.getShape()[0];
		
		minlat = latarr.getFloat(latindex.set(0));
		maxlat = latarr.getFloat(latindex.set(nlat-1));
		float rangelat = maxlat-minlat;
		
		minlon = lonarr.getFloat(lonindex.set(0));
		maxlon = lonarr.getFloat(lonindex.set(nlon-1));
		float rangelon = maxlon-minlon;
		
		// Keeping the same number of cells, what's the mean distance?
		
		float latdist = (rangelat)/nlat;
		float londist = (rangelon)/nlon;
		
		float cellsize = (latdist+londist)/2f;
		
		int new_nlat = (int) Math.floor(rangelat/cellsize);
		int new_nlon = (int) Math.floor(rangelon/cellsize);
		
		// Create new arrays for use in dimensions
		
		float[] lataxis = new float[new_nlat];
		float[] lonaxis = new float[new_nlon];
		
		float latc = minlat;
		float lonc = minlon;
	
		for(int j = 0; j < new_nlat;j++){
			lataxis[j] = latc;
			latc += cellsize;
		}
		
		for(int i = 0; i < new_nlon;i++){
			lonaxis[i] = lonc;
			lonc += cellsize;
		}
		
		// Ideally, if Variable has lat or lon as a dim, process, otherwise
		// write as is.
		
		writer = NetcdfFileWriter.createNew(NetcdfFileWriter.Version.netcdf4, outFile, null);
		
		List<Dimension> dims = new ArrayList<Dimension>();
		
		dims.add(writer.addDimension(null, outTimeName, tlen));
		dims.add(writer.addDimension(null, outLayerName, dlen));
		dims.add(writer.addDimension(null, outLatName, new_nlat));	
		dims.add(writer.addDimension(null, outLonName, new_nlon));
			
		
		outTime = writer.addVariable(null, outTimeName, DataType.FLOAT, dims
					.subList(0, 1));
		
		outLayer = writer.addVariable(null, outLayerName, DataType.FLOAT, dims
				.subList(1, 2));

		outLat = writer.addVariable(null, outLatName, DataType.FLOAT, dims
				.subList(2, 3));
		
		outLon = writer.addVariable(null, outLonName, DataType.FLOAT, dims
				.subList(3, 4));		
		
		outVar = writer.addVariable(null,outVarName, DataType.FLOAT, dims);
		
		writeAttributeInformation();
		
		writer.create();
		
		writer.write(outTime, inTime.read());
		writer.write(outLayer, inDepth.read());
		writer.write(outLat, Array.factory(float.class, new int[] {new_nlat}, lataxis));
		writer.write(outLon, Array.factory(float.class, new int[] {new_nlon}, lonaxis));
		
		float ilat,ilon;
		float ul,ll,ur,lr;
		double x,y,bs;
		int llat,llon,ulat,ulon;
		Array vararr;
		Array iarr;
		Index idx,odx;
		
		// Here, there would have to be some sorted of nested recursion for each variable...
		// and interpolate only along axes including lat and lon... ugh...
		
		for(int t = 0; t < tlen;t++){
			
			if(t%10 == 0){System.out.println("Processing step " + t + ". " + (int) ((double)t/(double)tlen*100) + "% complete.");}
			
			for(int k = 0; k < dlen;k++){
				vararr = inVar.read(new int[]{t,k,0,0}, new int[]{1,1,nlat,nlon});
				vararr = vararr.reduce();
				idx = vararr.getIndex();
				iarr = Array.factory(float.class, new int[] {1,1,new_nlat,new_nlon});
				odx = iarr.getIndex();
				
				for(int j = 0; j < new_nlat;j++){
					for(int i = 0; i < new_nlon;i++){
					
						// OK... first we need the latitude and longitude of the point to be interpolated...
						
						ilat = lataxis[j];
						ilon = lonaxis[i];
						
						llat = Arrays.binarySearch(platarr, ilat);
						if(llat<0){llat = -(llat+2);}
						ulat = llat+1;
						llon = Arrays.binarySearch(plonarr, ilon);
						if(llon<0){
							llon = -(llon+2);}
						ulon = llon+1;
						
						ul = vararr.getFloat(idx.set(ulat,llon));
						ur = vararr.getFloat(idx.set(ulat,ulon));
						ll = vararr.getFloat(idx.set(llat,llon));
						lr = vararr.getFloat(idx.set(llat,ulon));
						
						if(ul==NODATA||ur==NODATA||ll==NODATA||lr==NODATA){
							iarr.setFloat(odx.set(0,0,j,i), NODATA);
							
							// Wheeee!  Bilinear interpolation!
							
						}else{
							x = (ilon - plonarr[llon])/(plonarr[ulon]-plonarr[llon]);
							y = (ilat - platarr[llat])/(platarr[ulat]-platarr[llat]);
							bs = (ll*(1-x)*(1-y))+(lr*x*(1-y))+(ul*(1-x)*y)+(ur*x*y);
							iarr.setDouble(odx.set(0,0,j,i),bs);
						}					
					}
				}
				
				writer.write(outVar, new int[] {t,k,0,0}, iarr);
				
			}
		}
		
		writer.close();
	}
	
	/**
	 * Clones attribute information from source variables to destination variables.
	 * 
	 * @param infile - The input file
	 * @param inVarName - The name of the input variable
	 * @param outfile - The output file
	 * @param outVarName - The name of the output variable
	 */

	private void cloneAttributes(NetcdfFile infile, String inVarName,
			NetcdfFileWriter writer, Variable outVar) {

		// Find the variable
		
		Variable vi = infile.findVariable(inVarName);

		// Grab all of its attributes - unchecked, but should be OK.
		
		List<Attribute> l = vi.getAttributes();
		for (Attribute a : l) {

			writer.addVariableAttribute(outVar, a);
		}
	}
	
	private void writeAttributeInformation() {

		cloneAttributes(ncf, inTimeName, writer, outTime);
		cloneAttributes(ncf, inLayerName, writer, outLayer);
		writer.addVariableAttribute(outLon, new Attribute("units", "degrees_east"));
		writer.addVariableAttribute(outLon, new Attribute("standard_name", outLonName));
		writer.addVariableAttribute(outLon, new Attribute("axis", "x"));
		writer.addVariableAttribute(outLon, new Attribute("valid_range", Array.factory(
				float.class, new int[] { 2 }, new float[] { minlon, maxlon})));

		writer.addVariableAttribute(outLat, new Attribute("units", "degrees_north"));
		writer.addVariableAttribute(outLat, new Attribute("standard_name", outLatName));
		writer.addVariableAttribute(outLat, new Attribute("axis", "y"));
		writer.addVariableAttribute(outLat, new Attribute("valid_range", Array.factory(
				float.class, new int[] { 2 }, new float[] { minlat, maxlat })));
		writer.addVariableAttribute(outVar, new Attribute("units", "ms-1"));
		writer.addVariableAttribute(outVar, new Attribute("missing_value", 1E34));
		writer.addVariableAttribute(outVar, new Attribute("valid_range", Array.factory(
				float.class, new int[] { 2 },
				new float[] { lb_valid, ub_valid })));
	}
		
}
