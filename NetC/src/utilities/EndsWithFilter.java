package utilities;

import java.io.File;
import java.io.FilenameFilter;

/*
 * Helper Filter class - used for filtering out file names having the provided
 * extension.
 */

class EndsWithFilter implements FilenameFilter {
	private String extension;

	public EndsWithFilter(String extension) {
		this.extension = extension;
	}

	public boolean accept(File dir, String name) {
		if (name.endsWith(extension))
			return true;
		else
			return (new File(dir, name)).isDirectory();
	}
}
