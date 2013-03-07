package test;

import static org.junit.Assert.*;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import conv.MakeW;

public class TestMakeW {

	MakeW mw = new MakeW();
	
	public static void main(String[] args){
		TestMakeW tmw = new TestMakeW();
		tmw.setUp();
		tmw.testGo();
		System.out.println("Complete");
		
	}
	
	@Before
	public void setUp(){
		mw.setReproject(false);
		mw.setInputUFile("C:/Temp/Linear_X.nc");
		mw.setInputVFile("C:/Temp/Linear_X.nc");
		mw.setOutputWFile("C:/Temp/Test_W.nc");
		mw.setInputBathyFile("E:/HPC/Modeling/AUS/Input/NetCDF/GLB_depth_2d.nc");
		mw.setInTimeName("Time");
		mw.setOutTimeName("Time");
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testGo() {
		mw.go();
	}
}
