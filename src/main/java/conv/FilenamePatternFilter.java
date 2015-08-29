package conv;

import java.io.File;
import java.io.FilenameFilter;

/**
 * Helper Filter class - used for filtering out file names having the provided
 * extension.
 */

public class FilenamePatternFilter implements FilenameFilter {
	private String pattern;

	/**
	 * Constructor accepting string input representing the file extension to include.
	 * 
	 * @param extension
	 */
	
	public FilenamePatternFilter(String pattern) {
		this.pattern = pattern;
	}

	/**
	 * identifies whether the provided file has the appropriate extension
	 * 
	 * @param dir - The directory to be searched
	 * @param name - The name of the file to be queried
	 */
	
	public boolean accept(File dir, String name) {
		if (name.matches(pattern)){
			return true;}
		return false;
	}
}