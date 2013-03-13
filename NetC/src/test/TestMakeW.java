package test;

//import static org.junit.Assert.*;

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
		//mw.setReproject(false);
		//mw.setInputUFile("C:/Temp/Linear_Xn2_u.nc");
		//mw.setInputVFile("C:/Temp/Zeros_n2_v.nc");
		//mw.setInputBathyFile("C:/Temp/Floor2_Negative_1.nc");
		//mw.setOutputWFile("C:/Temp/Test_W.nc");
		//mw.setInTimeName("Time");
		//mw.setOutTimeName("Time");
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testGo() {
		mw.go();
	}
}
