package conv;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.util.StringTokenizer;

public class ASCIIGRID {

	private int ncols;
	private int nrows;
	private float xllcorner;
	private float yllcorner;
	private float cellsize;
	private float NODATA = -9999;
	private ByteBuffer vertices;
	private ByteBuffer indices;
	private float[] xdim;
	private float[] ydim;
	private int triangles;

	public ASCIIGRID() {
	}

	public ASCIIGRID(String path) {
		readFloatFile(path);
	}

	public void readFloatFile(String path) {
		BufferedReader br = null;

		// Read the header information

		try {
			br = new BufferedReader(new FileReader(path));
			StringTokenizer stk;
			String ln = br.readLine();
			stk = new StringTokenizer(ln);
			stk.nextToken();
			ncols = Integer.parseInt(stk.nextToken());
			ln = br.readLine();
			stk = new StringTokenizer(ln);
			stk.nextToken();
			nrows = Integer.parseInt(stk.nextToken());
			ln = br.readLine();
			stk = new StringTokenizer(ln);
			stk.nextToken();
			xllcorner = Float.parseFloat(stk.nextToken());
			ln = br.readLine();
			stk = new StringTokenizer(ln);
			stk.nextToken();
			yllcorner = Float.parseFloat(stk.nextToken());
			ln = br.readLine();
			stk = new StringTokenizer(ln);
			stk.nextToken();
			cellsize = Float.parseFloat(stk.nextToken());
			ln = br.readLine();
			stk = new StringTokenizer(ln);
			stk.nextToken();
			NODATA = Integer.parseInt(stk.nextToken());

			xdim = new float[ncols];
			ydim = new float[nrows];

			for (int i = 0; i < ncols; i++) {
				xdim[i] = xllcorner + (cellsize * (float) i);
			}

			for (int j = 0; j < nrows; j++) {
				ydim[j] = yllcorner + (cellsize * (float) j);
			}

			triangles = 2 * (nrows - 1) * (ncols - 1);
			
			vertices = ByteBuffer.allocateDirect(ncols * nrows * 3 * 4).order(ByteOrder.nativeOrder());
			indices = ByteBuffer.allocateDirect(triangles * 3 * 4).order(ByteOrder.nativeOrder());
			
			for (int i = nrows - 1; i >= 0; i--) {
				ln = br.readLine();
				stk = new StringTokenizer(ln);
				for (int j = 0; j < ncols; j++) {
					String s = stk.nextToken();
					float f = Float.parseFloat(s);
					int index = i + j * ncols;
					vertices.putFloat((index*3 + 0) * 4, xllcorner + ((float)j*cellsize));
					vertices.putFloat((index*3 + 1) * 4, yllcorner + ((float)j*cellsize));
					vertices.putFloat((index*3 + 2) * 4, f == NODATA ? Float.NaN : f);
				}
			}
			
			indices.clear();
			for (int i = 0; i < ncols - 1; i++) {
				for (int j = 0; j < nrows - 1; j++) {
					indices.putInt(j * ncols + i);
					indices.putInt(j * ncols + i + 1);
					indices.putInt((j + 1) * ncols + i + 1);

					indices.putInt(j * ncols + i);
					indices.putInt((j + 1) * ncols + i + 1);
					indices.putInt((j + 1) * ncols + i);
				}
			}
			
			indices.flip();

		} catch (NumberFormatException e) {
			e.printStackTrace();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		} finally {
			if (br != null) {
				try {
					br.close();
				} catch (IOException e) {
				}
			}
		}
	}

	float getFloatVertex(int i, int j) {
		return vertices.getFloat(((i * nrows) + j) * 4);
	}

	public int getNcols() {
		return ncols;
	}

	public void setNcols(int ncols) {
		this.ncols = ncols;
	}

	public int getNrows() {
		return nrows;
	}

	public void setNrows(int nrows) {
		this.nrows = nrows;
	}

	public float getXllcorner() {
		return xllcorner;
	}

	public void setXllcorner(float xllcorner) {
		this.xllcorner = xllcorner;
	}

	public float getYllcorner() {
		return yllcorner;
	}

	public void setYllcorner(float yllcorner) {
		this.yllcorner = yllcorner;
	}

	public float getCellsize() {
		return cellsize;
	}

	public void setCellsize(float cellsize) {
		this.cellsize = cellsize;
	}

	public ByteBuffer getVertices() {
		return vertices;
	}

	public void setVertices(ByteBuffer vertices) {
		this.vertices = vertices;
	}

	public ByteBuffer getIndices() {
		return indices;
	}
	
	public int getNumTriangles(){
		return triangles;
	}

	public void setIndices(ByteBuffer indices) {
		this.indices = indices;
	}
	
	public static void main(String[] args){
		ASCIIGRID asc = new ASCIIGRID("C:/Temp/xmap.txt");
		System.out.println();
		System.out.println(asc.getFloatVertex(0, 0));
		System.out.println(asc.getFloatVertex(0, 1));
		System.out.println(asc.getFloatVertex(1, 0));
		System.out.println(asc.getFloatVertex(0, 0));	
	}
}
