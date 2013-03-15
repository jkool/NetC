package test;

import org.junit.After;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import ucar.ma2.Array;
import ucar.ma2.Index;

import conv.MakeW;

public class TestMakeW {

	double e = 1E-16;

	MakeW mw = new MakeW();

	Array zeros = Array.factory(Double.class, new int[] { 3, 3 }, new double[] {
			0, 0, 0, 0, 0, 0, 0, 0, 0 });
	Array ones = Array.factory(Double.class, new int[] { 3, 3 }, new double[] {
			1, 1, 1, 1, 1, 1, 1, 1, 1 });
	Array u1 = Array.factory(Double.class, new int[] { 3, 3 }, new double[] {
			1, 2, 3, 1, 2, 3, 1, 2, 3 });
	Array v1 = Array.factory(Double.class, new int[] { 3, 3 }, new double[] {
			1, 1, 1, 2, 2, 2, 3, 3, 3 });
	Array u2 = Array.factory(Double.class, new int[] { 3, 3 }, new double[] {
			1, 3, 5, 1, 3, 5, 1, 3, 5 });
	Array v2 = Array.factory(Double.class, new int[] { 3, 3 }, new double[] {
			1, 1, 1, 3, 3, 3, 5, 5, 5 });
	Array u3 = Array.factory(Double.class, new int[] { 3, 3 }, new double[] {
			3, 2, 1, 3, 2, 1, 3, 2, 1 });
	Array v3 = Array.factory(Double.class, new int[] { 3, 3 }, new double[] {
			3, 3, 3, 2, 2, 2, 1, 1, 1 });

	Array zeros_3D;
	Array ones_3D;
	Array u3D1;
	Array v3D1;
	Array u3D2;
	Array v3D2;

	Array seq = Array.factory(Double.class, new int[] { 3, 3 }, new double[] {
			2, 4, 6, 12, 15, 18, 28, 32, 36 });

	Array lats1 = u1;
	Array lons1 = v1;
	Array lats2 = u2;
	Array lons2 = v2;
	Array lats3 = Array.factory(Double.class, new int[] { 3, 3 }, new double[] {
			0, 10, 20, 0, 10, 20, 0, 10, 20 });
	Array lons3 = Array.factory(Double.class, new int[] { 3, 3 }, new double[] {
			0, 0, 0, 5, 5, 5, 10, 10, 10 });

	float[] zero_arr = new float[] { 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	float[] ones_arr = new float[] { 0,  1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
	float[] seq_arr1 = new float[] { 0, -10, -9, -8, -7, -6, -5, -4, -3, -2, -1 };
	float[] seq_arr2 = new float[] { 0, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1 };
	float[] seq_arr3 = new float[] { 0, -20, -18, -16, -14, -12, -10, -8, -6, -4, -2};

	public static void main(String[] args) {
		TestMakeW tmw = new TestMakeW();
		// tmw.testLinearX();
		// tmw.testLinearY();
		// tmw.testSlopeX();

		System.out.println("Complete");
	}

	@Before
	public void setUp() {
		int nlayers = 11;
		zeros_3D = Array.factory(Double.class, new int[] { 1, nlayers, 3, 3 });
		ones_3D = Array.factory(Double.class, new int[] { 1, nlayers, 3, 3 });
		u3D1 = Array.factory(Double.class, new int[] { 1, nlayers, 3, 3 });
		v3D1 = Array.factory(Double.class, new int[] { 1, nlayers, 3, 3 });
		u3D2 = Array.factory(Double.class, new int[] { 1, nlayers, 3, 3 });
		v3D2 = Array.factory(Double.class, new int[] { 1, nlayers, 3, 3 });

		Index idx = zeros_3D.getIndex();
		for (int k = 0; k < nlayers; k++) {
			for (int i = 0; i < 3; i++) {
				for (int j = 0; j < 3; j++) {
					zeros_3D.setDouble(idx.set(0, k, i, j), 0);
					ones_3D.setDouble(idx, 1);
					u3D1.setDouble(idx, j);
					v3D1.setDouble(idx, i);
					u3D2.setDouble(idx, 2-j);
					v3D2.setDouble(idx, 2-i);
				}
			}
		}
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testLinearX() {
		mw.setReproject(false);
		mw.setInputUFile("C:/Temp/Linear_Xn2_u.nc");
		mw.setInputVFile("C:/Temp/Zeros_n2_v.nc");
		mw.setOutputWFile("C:/Temp/Test_W_x.nc");
		mw.setInTimeName("Time");
		mw.setOutTimeName("Time");
		mw.go();
	}

	@Test
	public void testLinearY() {
		mw.setReproject(false);
		mw.setInputUFile("C:/Temp/Zeros_n2_u.nc");
		mw.setInputVFile("C:/Temp/LinearY_n2_v.nc");
		mw.setOutputWFile("C:/Temp/Test_W_y.nc");
		mw.setInTimeName("Time");
		mw.setOutTimeName("Time");
		mw.go();
	}

	@Test
	public void testDx() {
		Assert.assertEquals(0.0d, mw.dx(zeros), e);
		Assert.assertEquals(0.0d, mw.dx(ones), e);
		Assert.assertEquals(1.0d, mw.dx(u1), e);
		Assert.assertEquals(0.0d, mw.dx(v1), e);
		Assert.assertEquals(2.0d, mw.dx(u2), e);
		Assert.assertEquals(0.0d, mw.dx(v2), e);
		Assert.assertEquals(-1.0d, mw.dx(u3), e);
		Assert.assertEquals(0.0d, mw.dx(v3), e);
		Assert.assertEquals(3.0d, mw.dx(seq), e);
	}

	@Test
	public void testDy() {
		Assert.assertEquals(0.0d, mw.dy(zeros), e);
		Assert.assertEquals(0.0d, mw.dy(ones), e);
		Assert.assertEquals(0.0d, mw.dy(u1), e);
		Assert.assertEquals(1.0d, mw.dy(v1), e);
		Assert.assertEquals(0.0d, mw.dy(u2), e);
		Assert.assertEquals(2.0d, mw.dy(v2), e);
		Assert.assertEquals(0.0d, mw.dy(u3), e);
		Assert.assertEquals(-1.0d, mw.dy(v3), e);
		Assert.assertEquals(14d, mw.dy(seq), e);
	}

	@Test
	public void testCalcDwdz() {
		mw.setReproject(false);
		Assert.assertEquals(-1, mw.calcdwdz(lats1, lons1, u1, zeros), e);
		Assert.assertEquals(-1, mw.calcdwdz(lats1, lons1, zeros, v1), e);
		Assert.assertEquals(-1, mw.calcdwdz(lats1, lons1, u1, ones), e);
		Assert.assertEquals(-1, mw.calcdwdz(lats1, lons1, ones, v1), e);
		Assert.assertEquals(-2, mw.calcdwdz(lats1, lons1, u1, v1), e);
		Assert.assertEquals(-0.5, mw.calcdwdz(lats2, lons1, u1, zeros), e);
		Assert.assertEquals(-1.0, mw.calcdwdz(lats1, lons2, u1, zeros), e);
		Assert.assertEquals(-0.5, mw.calcdwdz(lats2, lons2, u1, zeros), e);
		Assert.assertEquals(-1.0, mw.calcdwdz(lats2, lons1, zeros, v1), e);
		Assert.assertEquals(-0.5, mw.calcdwdz(lats1, lons2, zeros, v1), e);
		Assert.assertEquals(-0.5, mw.calcdwdz(lats2, lons2, zeros, v1), e);
		Assert.assertEquals(-1.0, mw.calcdwdz(lats2, lons2, u1, v1), e);
		Assert.assertEquals(1.0, mw.calcdwdz(lats2, lons2, u3, v3), e);
	}

	@Test
	public void testIntegrate() {
		mw.setReproject(false);
		Assert.assertArrayEquals(zero_arr, 
				mw.integrate(zeros_3D, zeros_3D, lons1, lats1), (float) e);
		Assert.assertArrayEquals(
				zero_arr, mw.integrate(ones_3D, zeros_3D, lons1, lats1), (float) e);
		Assert.assertArrayEquals(
				zero_arr, mw.integrate(zeros_3D, ones_3D, lons1, lats1), (float) e);
		Assert.assertArrayEquals(
				zero_arr, mw.integrate(ones_3D, ones_3D, lons1, lats1), (float) e);
		Assert.assertArrayEquals(
				seq_arr1, mw.integrate(u3D1, zeros_3D, lons1, lats1), (float) e);
		Assert.assertArrayEquals(
				seq_arr1, mw.integrate(zeros_3D, v3D1, lons1, lats1), (float) e);
		Assert.assertArrayEquals(
				zero_arr, mw.integrate(v3D1, u3D1, lons1, lats1), (float) e);
		Assert.assertArrayEquals(
				seq_arr3, mw.integrate(u3D1, v3D1, lons1, lats1), (float) e);
		Assert.assertArrayEquals(
				seq_arr1, mw.integrate(u3D1, v3D1, lons2, lats2), (float) e);
	}
}
