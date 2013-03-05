package conv;

import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;

/**
 * Performs basic time conversions to and from milliseconds.
 * 
 * @author Johnathan Kool
 *
 */

public class TimeConvert {

	public static final long SECS_IN_DAY = 24 * 60 * 60;
	public static final long MILLIS_IN_DAY = 24 * 60 * 60 * 1000;
	public static final long HYCOM_OFFSET = -2177539200000l;
	
	/**
	 * Converts time measurements into milliseconds.
	 * 
	 * @param unit - Original units of time
	 * @param val - Number of units
	 * @return
	 */

	public static long convertToMillis(String unit, String val) {

		if (unit.equalsIgnoreCase("Days") || unit.equalsIgnoreCase("Day")
				|| unit.equalsIgnoreCase("D")) {
			return daysToMillis(val);
		} else if (unit.equalsIgnoreCase("Seconds")
				|| unit.equalsIgnoreCase("Second")
				|| unit.equalsIgnoreCase("Secs")
				|| unit.equalsIgnoreCase("Sec") 
				|| unit.equalsIgnoreCase("S")) {
			return secondsToMillis(val);
		} else if (unit.equalsIgnoreCase("Date")) {
			try {
				DateFormat df = new SimpleDateFormat("M/d/yyyy");
				Date d = df.parse(val);
				return d.getTime();
			} catch (ParseException e) {
				throw new IllegalArgumentException("Date could not be parsed: "
						+ val + ".  Required format is M/d/yyyy.");
			}
			} else if (unit.equalsIgnoreCase("Hours")
				|| unit.equalsIgnoreCase("Hour")
				|| unit.equalsIgnoreCase("Hrs") 
				|| unit.equalsIgnoreCase("Hr")
				|| unit.equalsIgnoreCase("H")) {
			return hoursToMillis(val);
		} else if (unit.equalsIgnoreCase("Minutes")
				|| unit.equalsIgnoreCase("Minute")
				|| unit.equalsIgnoreCase("Mins")
				|| unit.equalsIgnoreCase("Min") 
				|| unit.equalsIgnoreCase("M")) {
			return minutesToMillis(val);
		}

		else
			throw new UnsupportedOperationException("Conversion of unit "
					+ unit + "has not been implemented");

	}

	/**
	 * Converts time measurements to milliseconds.
	 * 
	 * @param unit - Original units of time
	 * @param val - Number of units (as a double)
	 * @return
	 */

	
	public static long convertToMillis(String unit, double val) {

		if (unit.equalsIgnoreCase("Days") || unit.equalsIgnoreCase("Day")
				|| unit.equalsIgnoreCase("D")) {
			return daysToMillis(val);
		} else if (unit.equalsIgnoreCase("Seconds")
				|| unit.equalsIgnoreCase("Second")
				|| unit.equalsIgnoreCase("Secs")
				|| unit.equalsIgnoreCase("Sec") || unit.equalsIgnoreCase("S")) {
			return secondsToMillis(val);
		} else if (unit.equalsIgnoreCase("Hours")
				|| unit.equalsIgnoreCase("Hour")
				|| unit.equalsIgnoreCase("Hrs") || unit.equalsIgnoreCase("Hr")
				|| unit.equalsIgnoreCase("H")) {
			return hoursToMillis(val);
		} else if (unit.equalsIgnoreCase("Minutes")
				|| unit.equalsIgnoreCase("Minute")
				|| unit.equalsIgnoreCase("Mins")
				|| unit.equalsIgnoreCase("Min") || unit.equalsIgnoreCase("M")) {
			return minutesToMillis(val);
		} else if (unit.equalsIgnoreCase("Milliseconds")
				|| unit.equalsIgnoreCase("Millisecond")
				|| unit.equalsIgnoreCase("ms")) {
			return (long) val;
		} else
			throw new UnsupportedOperationException("Conversion of unit "
					+ unit + "has not been implemented");

	}

	/**
	 * Converts days into milliseconds
	 * 
	 * @param days - Number of days
	 * @return
	 */
	
	public static long daysToMillis(String days) {
		double d = Double.parseDouble(days);
		return Math.round(d * 24 * 60 * 60 * 1000);
	}

	/**
	 * Convenience function for Converting seconds into milliseconds 
	 * (multiply by 1000).
	 * 
	 * @param seconds
	 * @return
	 */
	
	public static long secondsToMillis(String seconds) {
		double d = Double.parseDouble(seconds);
		return Math.round(d * 1000);
	}

	/**
	 * Convert a date into milliseconds using the specified format.
	 * 
	 * @param date - String representation of the date to be converted.
	 * @return
	 */
	
	public static long dateToMillis(String date) {
		DateFormat df = new SimpleDateFormat("M/d/yyyy");
		Date d;
		try {
			d = df.parse(date);
			return d.getTime();
		} catch (ParseException e) {
			System.out.println("WARNING:  Date provided: " + date
					+ "could not be parsed.\n\n");
			e.printStackTrace();
		}
		return -1;
	}
	
	/**
	 * Converts hours into milliseconds
	 * 
	 * @param hours - String representation of the number of hours to be converted.
	 * @return
	 */

	public static long hoursToMillis(String hours) {
		double d = Double.parseDouble(hours);
		return Math.round(d * 60 * 60 * 1000);
	}

	/**
	 * Converts minutes into milliseconds
	 * 
	 * @param minutes - String representation of the number of minutes to be converted.
	 * @return
	 */
	
	public static long minutesToMillis(String minutes) {
		double d = Double.parseDouble(minutes);
		return Math.round(d * 60 * 1000);
	}

	/**
	 * Converts days into milliseconds
	 * 
	 * @param days - String representation of the number of days to be converted.
	 * @return
	 */
	
	public static long daysToMillis(double days) {
		return Math.round(days * 24 * 60 * 60 * 1000);
	}

	/**
	 * Converts from seconds to milliseconds
	 * 
	 * @param seconds
	 * @return
	 */
	
	public static long secondsToMillis(double seconds) {
		return Math.round(seconds * 1000);
	}

	/**
	 * Converts from hours to milliseconds
	 * 
	 * @param hours
	 * @return
	 */
	
	public static long hoursToMillis(double hours) {
		return Math.round(hours * 60 * 60 * 1000);
	}

	/**
	 * Converts from minutes to milliseconds
	 * 
	 * @param minutes
	 * @return
	 */
	
	public static long minutesToMillis(double minutes) {
		return Math.round(minutes * 60 * 1000);
	}
	
	/**
	 * Converts a HYCOM time to Java milliseconds
	 * 
	 * @param days
	 * @return
	 */
	
	public static long HYCOMToMillis(long days){
		return daysToMillis(days)+HYCOM_OFFSET;
	}
	
	/**
	 * Converts Java milliseconds to a HYCOM time 
	 * @return
	 */
	
	public static double millisToHYCOM(double millis){
		return millisToDays(millis-HYCOM_OFFSET);
	}

	/** 
	 * Converts a Java long value (days) to a HYCOM time 
	 * @return
	 */
	
	public static long javaLongToHYCOM(long days){
		return days - (long) millisToDays(-HYCOM_OFFSET);
	}

	/** 
	 * Converts a Java Date to a HYCOM time 
	 * @return
	 */
	
	public static long javaDateToHYCOM(Date d){
		return (long) millisToDays(d.getTime()-HYCOM_OFFSET);
	}
	
	public static long HYCOMToJavaLong(long days){
		return days + (long) millisToDays(-HYCOM_OFFSET);
	}
	
	public static Date HYCOMToJavaDate(long days){
		return new Date(daysToMillis(days)+HYCOM_OFFSET);
	}

	public static double convertFromMillis(String unit, double millis) {
		if (unit.equalsIgnoreCase("Days") || unit.equalsIgnoreCase("Day")
				|| unit.equalsIgnoreCase("D")) {
			return millisToDays(millis);
		} else if (unit.equalsIgnoreCase("Seconds")
				|| unit.equalsIgnoreCase("Second")
				|| unit.equalsIgnoreCase("Secs")
				|| unit.equalsIgnoreCase("Sec") || unit.equalsIgnoreCase("S")) {
			return millisToSeconds(millis);
		} else if (unit.equalsIgnoreCase("Hours")
				|| unit.equalsIgnoreCase("Hour")
				|| unit.equalsIgnoreCase("Hrs") || unit.equalsIgnoreCase("Hr")
				|| unit.equalsIgnoreCase("H")) {
			return millisToHours(millis);
		} else if (unit.equalsIgnoreCase("Minutes")
				|| unit.equalsIgnoreCase("Minute")
				|| unit.equalsIgnoreCase("Mins")
				|| unit.equalsIgnoreCase("Min") || unit.equalsIgnoreCase("M")) {
			return millisToMinutes(millis);
		}

		else
			throw new UnsupportedOperationException("Conversion of unit "
					+ unit + "has not been implemented");

	}

	public static double millisToDays(double millis) {
		return (millis) / (24f * 60f * 60f * 1000f);
	}
	
	public static double millisToDays(long millis) {
		return ((double) millis) / (24f * 60f * 60f * 1000f);
	}

	public static double millisToSeconds(double millis) {
		return millis / 1000f;
	}

	public static double millisToHours(double millis) {
		return millis / (60f * 60f * 1000f);
	}

	public static double millisToMinutes(double millis) {
		return millis / (60f * 1000f);
	}
	
	public static Date millisToDate(double millis){
		return new Date((long) millis);
	}
	
	/**
	 * Converts from milliseconds to a simple String representation.
	 * 
	 * @param millis
	 * @return
	 */
	
	public static String millisToString(double millis){
		double remainder = millis;
		double day = 1000*60*60*24;
		int days = (int) Math.floor(millis/day);
		remainder = millis-(days*day);
		double hour = 1000*60*60;
		int hours = (int) Math.floor(remainder/hour);
		remainder = remainder-(hours*hour);
		double minute = 1000*60;
		int minutes = (int) Math.floor(remainder/minute);
		remainder = remainder-(minute*minutes);
		double seconds = remainder/1000;
		
		StringBuffer sb = new StringBuffer();

		if(days!=0){
			sb.append(days + "d, ");
		}
		if(hours!=0 || days!=0){
			sb.append(hours + "h ");
		}
		if(minutes!=0 || days!=0 ||hours!=0){
			sb.append(minutes + "m ");
		}
		sb.append(seconds + "s");

		return sb.toString();	
	}
}
