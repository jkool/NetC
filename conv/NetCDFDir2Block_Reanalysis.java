package conv;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.ParseException;
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

/**
 * Used to convert single NetCDF time slices into processed blocks.
 * 
 * @author Johnathan Kool
 */

public class NetCDFDir2Block_Reanalysis {

    private double maxFileSize = 4E9;
    // private String fileDir = "F:/MH";
    private String fileDir = "G:\\HI";
    private String inVarName = "u";
    private String inHorizontalDim = "lon";
    private String inVerticalDim = "lat";
    private String inTimeName = "time";
    private String inDepthName = "depth";
    private String inLatName = "lat";
    private String inLonName = "lon";
    private String outTimeName = "Time";
    private String outDepthName = "Depth";
    private String outLatName = "Latitude";
    private String outLonName = "Longitude";
    private String outVarName = inVarName;
    // private String outputDir = "F:/MH_Blocks";
    private String outputDir = "G:\\Temp\\HYCOM\\HI\\blk";
    private DateFormat df = new SimpleDateFormat("yyyy_MM_dd");
    private DateFormat df2 = new SimpleDateFormat("yyyy_MM_dd zzz");
    private FilenameFilter filter = new FilenamePatternFilter(".*_"
            + inVarName + "_.*\\.nc");
    private int blockdays = 60;
    private long blocksize = TimeConvert.daysToMillis(Integer
            .toString(blockdays));
    private long initialTime;
    private long endTime = Long.MAX_VALUE;
    private Array timeArr;
    private Array depthArr;
    private Array latArr;
    private Array lonArr;
    private boolean reverseDepth = true;

    public static void main(String[] args) {
        NetCDFDir2Block_Reanalysis ncb = new NetCDFDir2Block_Reanalysis();
        
        if(args.length!=5){
            System.out.println("Usage: <input variable name> <input directory> <output directory> <block size> <start time yyyy_MM_dd in UTC>");
            System.exit(-1);
        }
        
        ncb.inVarName = args[0];
        ncb.fileDir = args[1];
        ncb.outputDir = args[2];
        ncb.blocksize= TimeConvert.daysToMillis(args[3]);
        ncb.filter = new FilenamePatternFilter(".*_"
                + ncb.inVarName + "_.*\\.nc");
        if (ncb.inVarName.equalsIgnoreCase("w")){
            ncb.inVarName="w_velocity";
        }
        if (ncb.inVarName.equalsIgnoreCase("u")){
            ncb.inVarName="water_u";
        }
        if (ncb.inVarName.equalsIgnoreCase("v")){
            ncb.inVarName="water_v";
        }
        String initString = args[4] + " UTC";
        //String initString = "1993_01_01" + " UTC";

        try {
            // ncb.initialTime = ncb.df2.parse("2014_03_06 UTC").getTime();
            ncb.initialTime = ncb.df2.parse(initString).getTime();
        } catch (ParseException e) {
            e.printStackTrace();
        }
        ncb.convert(ncb.fileDir, ncb.filter, ncb.blocksize);
        System.out.println("Complete.");
    }

    private NetcdfFile ncf_in;

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

    /**
     * Clones the attributes of one NetCDF file into another.
     * 
     * @param infile
     *            - The source NetCDF file
     * @param inVarName
     *            - The variable name to be copied
     * @param outfile
     *            - The destination NetCDF file
     * @param outVarName
     *            - The output variable name
     */

    private void cloneAttributes(NetcdfFile infile, String inVarName,
            NetcdfFileWriteable outfile, String outVarName) {

        // Find the variable

        Variable vi = infile.findVariable(inVarName);

        // Grab all of its attributes - unchecked, but should be OK.

        List<Attribute> l = vi.getAttributes();
        for (Attribute a : l) {
        	if(a.getName().equalsIgnoreCase("units") && inVarName.equalsIgnoreCase("time")){
        		Attribute aa = new Attribute("units","days since 1900-12-31 00:00:00");
        		outfile.addVariableAttribute(outVarName,aa);
        	}
        	
        	else if(a.getName().equalsIgnoreCase("time_origin") && inVarName.equalsIgnoreCase("time")){
        		Attribute aa = new Attribute("time_origin", "1900-12-31 00:00:00");
        		outfile.addVariableAttribute(outVarName,aa);
        	}
        	
        	else if(a.getName().equalsIgnoreCase("calendar")){
        		Attribute aa = new Attribute("calendar", "standard");
        		outfile.addVariableAttribute(outVarName,aa);        		
        	}
        	
        	else{
        		outfile.addVariableAttribute(outVarName, a);
        	}
        }
    }

    /**
     * Closes the output file
     * 
     * @param ncf_out
     */

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

    /**
     * Performs conversion of a directory of files into blocks
     */

    public void convert() {
        NetCDFDir ncdir = new NetCDFDir(fileDir, filter, inTimeName);
        ArrayList<Long> times = new ArrayList<Long>(ncdir.getTimes());

        Iterator<Long> it = times.iterator();

        while (it.hasNext()) {
            long t = it.next();
            if (t < initialTime || t > endTime) {
                it.remove();
            }
        }
        
        if(times.isEmpty()){System.out.println("WARNING: Number of filtered files is 0.");return;}

        List<Long> blocks = calcBlock(times, blocksize);
        // if (blocks.size()!=blockdays){
        // System.out.println("WARNING:  Block size is different than expected: "
        // + blocks.size() + " found versus " + blockdays + " expected.");
        // }
        TreeMap<Long, NetcdfFile> selection;

        // iterate through the blocks
        for (int i = 1; i < blocks.size(); i++) {
            System.out.println("Working on block " + i + "...");
            long start_time = blocks.get(i - 1);
            long end_time = blocks.get(i);
            selection = selectBetween(start_time, end_time, ncdir.getFiles());
            String varout = inVarName.equalsIgnoreCase("w_velocity")?"w":inVarName;
            varout = inVarName.equalsIgnoreCase("water_u")?"u":inVarName;
            varout = inVarName.equalsIgnoreCase("water_v")?"v":inVarName;
            outVarName = varout;
            String outputPath = outputDir + "/"
                    + df.format(new Date(start_time)) + "_to_"
                    + df.format(new Date(end_time)) + "_" + varout + ".nc";
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

    /**
     * Conversion method accepting function parameters
     * 
     * @param fileDir
     *            - the directory containing the files to be processed
     * @param filter
     *            - used to filter set of files for processing
     * @param blocksize
     *            - sets the size of the block
     */

    public void convert(String fileDir, FilenameFilter filter, long blocksize) {
        this.fileDir = fileDir;
        this.filter = filter;
        this.blocksize = blocksize;
        convert();
    }

    /**
     * Copies the contents of single NetCDF files into a block file
     * 
     * @param map
     * @param ncf_out
     */

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
            e.printStackTrace();
        } catch (InvalidRangeException e) {
            e.printStackTrace();
        }
    }

    /**
     * Generates an empty template for a block
     * 
     * @param files
     * @param outputFile
     * @return
     */

    private NetcdfFileWriteable makeTemplate(TreeMap<Long, NetcdfFile> files,
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
            int lonidx = lon.findDimensionIndex(inHorizontalDim);
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
                if (reverseDepth) {
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
            e.printStackTrace();
        } catch (InvalidRangeException e) {
            e.printStackTrace();
        }

        return ncf_out;
    }

    /**
     * Creates an array of ones of arbitrary length
     * 
     * @param length
     * @return
     */

    private int[] ones(int length) {
        int[] out = new int[length];
        for (int i = 0; i < length; i++) {
            out[i] = 1;
        }
        return out;
    }

    /**
     * Selects netCDF files after a given point in time (in Java milliseconds)
     * 
     * @param from
     *            - the time in milliseconds
     * @param map
     *            - the map from which files will be selected with associated
     *            time values as keys.
     * @return
     */

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

    /**
     * Selects a block of time from the entire list of files
     * 
     * @param from
     *            - Starting time (as Java milliseconds)
     * @param until
     *            - Ending time (up to, not including as java milliseconds)
     * @param map
     *            - The Map object containing the complete file list linked with
     *            time
     * @return
     */

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

    /**
     * Gets the size of the block
     * 
     * @return
     */

    public long getBlocksize() {
        return blocksize;
    }

    /**
     * Retrieves the path of the file directory for processing
     * 
     * @return
     */

    public String getFileDir() {
        return fileDir;
    }

    /**
     * Retrieves the FilenameFilter object being used to filter files.
     * 
     * @return
     */

    public FilenameFilter getFilter() {
        return filter;
    }

    /**
     * Retrieves the name of the input Depth variable
     * 
     * @return
     */

    public String getInDepthName() {
        return inDepthName;
    }

    /**
     * Retrieves the name of the input horizontal dimension (e.g. X)
     * 
     * @return
     */

    public String getInHorizontalDim() {
        return inHorizontalDim;
    }

    /**
     * Retrieves the initial time (start time for blocks)
     * 
     * @return
     */

    public long getInitialTime() {
        return initialTime;
    }

    /**
     * Retrieves the initial time in the provided units (start time for blocks)
     * 
     * @return
     */

    public double getInitialTime(String units) {
        return TimeConvert.convertFromMillis(units, initialTime);
    }

    /**
     * Retrieves the end time (start time for blocks)
     * 
     * @return
     */

    public long getEndTime() {
        return endTime;
    }

    /**
     * Retrieves the end time in the provided units (start time for blocks)
     * 
     * @return
     */

    public double getEndTime(String units) {
        return TimeConvert.convertFromMillis(units, initialTime);
    }

    /**
     * Retrieves the name of the input latitude variable
     * 
     * @return
     */

    public String getInLatName() {
        return inLatName;
    }

    /**
     * Retrieves the name of the input longitude variable
     * 
     * @return
     */

    public String getInLonName() {
        return inLonName;
    }

    /**
     * Retrieves the name of the input time variable
     * 
     * @return
     */

    public String getInTimeName() {
        return inTimeName;
    }

    /**
     * Retrieves the name of the input variable (e.g. u)
     * 
     * @return
     */

    public String getInVarName() {
        return inVarName;
    }

    /**
     * Retrieves the name of the input vertical dimension (e.g. Y)
     * 
     * @return
     */

    public String getInVerticalDim() {
        return inVerticalDim;
    }

    /**
     * Retrieves the maximum block file size (unused)
     * 
     * @return
     */

    public double getMaxFileSize() {
        return maxFileSize;
    }

    /**
     * Retrieves the name of the output Depth variable
     * 
     * @return
     */

    public String getOutDepthName() {
        return outDepthName;
    }

    /**
     * Retrieves the name of the output latitude variable
     * 
     * @return
     */

    public String getOutLatName() {
        return outLatName;
    }

    /**
     * Retrieves the name of the output longitude variable
     * 
     * @return
     */

    public String getOutLonName() {
        return outLonName;
    }

    /**
     * Retrieves the name of the directory that will contain the output block
     * files
     * 
     * @return
     */

    public String getOutputDir() {
        return outputDir;
    }

    /**
     * Retrieves the name of the output time variable
     * 
     * @return
     */

    public String getOutTimeName() {
        return outTimeName;
    }

    /**
     * Retrieves the name for the output variable
     * 
     * @return
     */

    public String getOutVarName() {
        return outVarName;
    }

    /**
     * Identifies whether depth values should be reversed (i.e. convert positive
     * values to negative)
     * 
     * @return
     */

    public boolean hasReverseDepth() {
        return reverseDepth;
    }

    /**
     * Sets the number of files to be aggregated into a block
     * 
     * @param blocksize
     */

    public void setBlocksize(long blocksize) {
        this.blocksize = blocksize;
    }

    /**
     * Sets the path of the file directory for processing
     * 
     * @param fileDir
     */

    public void setFileDir(String fileDir) {
        this.fileDir = fileDir;
    }

    /**
     * Sets the FilenameFilter object being used to filter files.
     * 
     * @param filter
     */

    public void setFilter(FilenameFilter filter) {
        this.filter = filter;
    }

    /**
     * Sets the name of the input Depth variable
     * 
     * @param inDepthName
     */

    public void setInDepthName(String inDepthName) {
        this.inDepthName = inDepthName;
    }

    /**
     * Sets the name of the input horizontal dimension (e.g. X)
     * 
     * @param inDepthName
     */

    public void setInHorizontalDim(String inHorizontalDim) {
        this.inHorizontalDim = inHorizontalDim;
    }

    /**
     * Sets the initial time for processing blocks
     * 
     * @param time
     */

    public void setInitialTime(long time) {
        this.initialTime = time;
    }

    /**
     * Sets the initial time for processing blocks using the specified units
     * 
     * @param time
     */

    public void setInitialTime(double val, String unit) {
        this.initialTime = TimeConvert.convertToMillis(unit, val);
    }

    /**
     * Sets the end time for processing blocks
     * 
     * @param time
     */

    public void setEndTime(long time) {
        this.endTime = time;
    }

    /**
     * Sets the end time for processing blocks using the specified units
     * 
     * @param time
     */

    public void setEndTime(double val, String unit) {
        this.endTime = TimeConvert.convertToMillis(unit, val);
    }

    /**
     * Sets the name of the input latitude variable
     * 
     * @param inLatName
     */

    public void setInLatName(String inLatName) {
        this.inLatName = inLatName;
    }

    /**
     * Sets the name of the input longitude variable
     * 
     * @param inLatName
     */

    public void setInLonName(String inLonName) {
        this.inLonName = inLonName;
    }

    /**
     * Sets the name of the input time variable
     * 
     * @param inLatName
     */

    public void setInTimeName(String inTimeName) {
        this.inTimeName = inTimeName;
    }

    /**
     * Sets the name of the input variable (e.g. u)
     * 
     * @param inLatName
     */

    public void setInVarName(String inVarName) {
        this.inVarName = inVarName;
    }

    /**
     * Sets the name of the input vertical dimension (e.g. Y)
     * 
     * @param inLatName
     */

    public void setInVerticalDim(String inVerticalDim) {
        this.inVerticalDim = inVerticalDim;
    }

    /**
     * Sets the maximum file size of the output block (unused)
     * 
     * @param maxFileSize
     */

    public void setMaxFileSize(double maxFileSize) {
        this.maxFileSize = maxFileSize;
    }

    /**
     * Sets whether depth values should be reversed (i.e. convert positive
     * values to negative)
     * 
     * @param negativeDepth
     */

    public void setReverseDepth(boolean negativeDepth) {
        this.reverseDepth = negativeDepth;
    }

    /**
     * Sets the name of the output depth variable
     * 
     * @param outDepthName
     */

    public void setOutDepthName(String outDepthName) {
        this.outDepthName = outDepthName;
    }

    /**
     * Sets the name of the output latitude variable
     * 
     * @param outLatName
     */

    public void setOutLatName(String outLatName) {
        this.outLatName = outLatName;
    }

    /**
     * Sets the name of the output longitude variable
     * 
     * @param outLonName
     */

    public void setOutLonName(String outLonName) {
        this.outLonName = outLonName;
    }

    /**
     * Sets the name of the directory to contain the output block files
     * 
     * @param outLatName
     */

    public void setOutputDir(String outputDir) {
        this.outputDir = outputDir;
    }

    /**
     * Sets the name of the output time variable
     * 
     * @param outLatName
     */

    public void setOutTimeName(String outTimeName) {
        this.outTimeName = outTimeName;
    }

    /**
     * Sets the name of the output variable (e.g. u)
     * 
     * @param outLatName
     */

    public void setOutVarName(String outVarName) {
        this.outVarName = outVarName;
    }
}

