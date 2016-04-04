#pragma once
#include "comm.h"

namespace biostacs_bmp {
#define __bit_depth 24

	static const char __MagicBitmapNum[11][8] =  //FONTSIZE 8
	{
		{ 0x00,0x18,0x24,0x24,0x24,0x24,0x24,0x18 }, //0
		{ 0x00,0x18,0x1c,0x18,0x18,0x18,0x18,0x18 }, //1
		{ 0x00,0x1e,0x30,0x30,0x1c,0x06,0x06,0x3e }, //2
		{ 0x00,0x1e,0x30,0x30,0x1c,0x30,0x30,0x1e }, //3
		{ 0x00,0x30,0x38,0x34,0x32,0x3e,0x30,0x30 }, //4
		{ 0x00,0x1e,0x02,0x1e,0x30,0x30,0x30,0x1e }, //5
		{ 0x00,0x1c,0x06,0x1e,0x36,0x36,0x36,0x1c }, //6
		{ 0x00,0x3f,0x30,0x18,0x18,0x0c,0x0c,0x0c }, //7
		{ 0x00,0x1c,0x36,0x36,0x1c,0x36,0x36,0x1c }, //8
		{ 0x00,0x1c,0x36,0x36,0x36,0x3c,0x30,0x1c }, //9
		{ 0x00,0x00,0x00,0x00,0x01,0x00,0x00,0x00 }//.
	};


	class BITMATWRITE500 {
	private:
		std::vector<std::vector<int> > datamatrix;
		std::string filename;
		std::vector<int> xval;
		std::vector<int> yval;
		inline int rowchk(int i) { return (((int)((i*__bit_depth + 31) / 32)) * 4); }
		inline int aabs(int x) { if (x<0) return (-x);  return(x); }
		inline int shrink(int i, int lower, double ratio) { double re = lower + double(i)*ratio;  return floor(re); }
		void NumOnText(int xt, int yt, char text);
		void line(int x0, int y0, int x1, int y1, int col);
		void drawframe();
		void shrinkline(int x0, int y0, int x1, int y1, int col);
		void PlotStepFunction(std::vector<int>x, std::vector<int>y);
		int widethoffset;
	public:
		BITMATWRITE500(
			std::string filename_,
			std::vector<int> xval_,
			std::vector<int> yval_) {
			std::vector<std::vector<int> > datamatrix2(500, std::vector<int>(500));
			datamatrix = datamatrix2;
			filename = filename_;
			xval = xval_;
			yval = yval_;
			widethoffset = 8 * 3;
		}
		void BITMAPCREATING_CT();

	};



	void BITMATWRITE500::NumOnText(int xt, int yt, char text) {
		int num = text - '0';
		int i, j, x;
		int b[8];
		if (text == '.') num = 10;

		for (i = 0; i<8; i++)
		{
			x = __MagicBitmapNum[num][i];
			for (j = 0; j<8; j++) {
				b[j] = x % 2;
				x = x / 2;
			}
			for (j = 0; j<8; j++) {
				if (b[j] % 2)
					datamatrix[yt + 8 - i][xt + j] = 1;
			}
		}
	}



	void BITMATWRITE500::line(int x0, int y0, int x1, int y1, int col) {
		//int n=(int)datamatrix[0].size();
		int dx = aabs(x1 - x0), sx = x0<x1 ? 1 : -1;
		int dy = aabs(y1 - y0), sy = y0<y1 ? 1 : -1;
		int err = (dx>dy ? dx : -dy) / 2, e2;
		for (;;) {
			datamatrix[y0][x0] = col;
			if (x0 == x1 && y0 == y1) break;
			e2 = err;
			if (e2 >-dx) { err -= dy; x0 += sx; }
			if (e2 < dy) { err += dx; y0 += sy; }
		}
	}

	void BITMATWRITE500::drawframe() {
		line(widethoffset, widethoffset, widethoffset, 500 - widethoffset, 1);
		line(widethoffset, 500 - widethoffset, 500 - widethoffset, 500 - widethoffset, 1);
		line(500 - widethoffset, 500 - widethoffset, 500 - widethoffset, widethoffset, 1);
		line(500 - widethoffset, widethoffset, widethoffset, widethoffset, 1);

		line(widethoffset, 250, widethoffset + 3, 250, 1);
		line(250, widethoffset, 250, widethoffset + 3, 1);

		/*Y X 0.5*/
		NumOnText(2, 245, '0');
		NumOnText(8, 245, '.');
		NumOnText(14, 245, '5');
		NumOnText(242, 2, '0');
		NumOnText(251, 2, '.');
		NumOnText(257, 2, '5');

		/*X Y0.25*/
		int mid = floor((widethoffset + 250) / 2);
		line(mid, widethoffset, mid, widethoffset + 2, 1);
		NumOnText(mid - 9, 2, '0');
		NumOnText(mid - 3, 2, '.');
		NumOnText(mid + 3, 2, '2');
		NumOnText(mid + 9, 2, '5');

		line(widethoffset, mid, widethoffset + 2, mid, 1);
		NumOnText(2, mid, '0');
		NumOnText(8, mid, '.');
		NumOnText(14, mid, '2');
		NumOnText(18, mid, '5');


		/*X 0.75*/
		mid = floor((-widethoffset + 250 + 500) / 2);
		line(mid, widethoffset, mid, widethoffset + 2, 1);
		NumOnText(mid - 9, 2, '0');
		NumOnText(mid - 3, 2, '.');
		NumOnText(mid + 3, 2, '7');
		NumOnText(mid + 9, 2, '5');
		line(widethoffset, mid, widethoffset + 2, mid, 1);
		NumOnText(2, mid, '0');
		NumOnText(8, mid, '.');
		NumOnText(14, mid, '7');
		NumOnText(18, mid, '5');

		/*X Y 1.0*/
		mid = 500 - widethoffset;
		line(mid, widethoffset, mid, widethoffset + 2, 1);
		NumOnText(mid - 3, 2, '1');
		NumOnText(mid + 3, 2, '.');
		NumOnText(mid + 9, 2, '0');

		NumOnText(2, mid, '1');
		NumOnText(8, mid, '.');
		NumOnText(14, mid, '0');
		line(mid, widethoffset, mid, widethoffset + 2, 1);

		/*origin: 0,0*/
		mid = widethoffset;

		NumOnText(mid - 3, 2, '0');
		NumOnText(mid + 3, 2, '.');
		NumOnText(mid + 9, 2, '0');



	}




	void BITMATWRITE500::shrinkline(int x0, int y0, int x1, int y1, int col) {
		double ratio = 1 - 2.*widethoffset / 500.;
		x0 = shrink(x0, widethoffset, ratio);
		y0 = shrink(y0, widethoffset, ratio);
		x1 = shrink(x1, widethoffset, ratio);
		y1 = shrink(y1, widethoffset, ratio);
		line(x0, y0, x1, y1, col);
	}


	void BITMATWRITE500::PlotStepFunction(std::vector<int>x, std::vector<int>y) {
		int n = x.size();
		shrinkline(1, 1, x[0], y[0], 2);
		for (int i = 1; i< n; i++) {
			shrinkline(x[i - 1], y[i - 1], x[i], y[i], 2);
		}
	}

	void BITMATWRITE500::BITMAPCREATING_CT() {
		const unsigned int offset = 54;
		std::ofstream fs;
		fs.open(filename.c_str(), std::ios::trunc);

		/******Diag********/
		shrinkline(1, 1, 500, 500, 3);

		/*CURVE*/
		PlotStepFunction(xval, yval);

		/**FRAME**/
		drawframe();

		unsigned char HEADR[offset + 4] = {
			0x42,0x4d,0xe6,0x71,0x0b,0x00,0x00,/*BW BMP*/
			0x00,0x00,0x00,
			0x36,0x00,0x00,0x00,
			0x28,0x00,0x00,0x00,
			0xf4,0x01,0x00,0x00,/*width: 500*/
			0xf4,0x01,0x00,0x00,/*heigth: 500*/
			0x01,0x00,0x18,
			0x00,0x00,0x00,0x00,0x00,
			0xb0,0x71,
			0x0b,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
			0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,
			0x00,0x00
		};
		//COLOR
		const unsigned char BLACK[3] = { 0x00,0x00,0x00 };
		const unsigned char WHITE[3] = { 0xff,0xff,0xff };
		const unsigned char RED[3] = { 0x00,0x00,0xff };
		const unsigned char GREEN[3] = { 0xcc,0xff,0xcc };

		//  unsigned char widthoff[1]={0x00};
		for (int i = 0; i<offset; i++) fs << HEADR[i];
		for (int i = 0; i<500; i++) {
			for (int j = 0; j<500; j++) {
				switch (datamatrix[i][j]) {
				case(1) :
					fs << BLACK[0] << BLACK[1] << BLACK[2];
					break;
				case(2) :
					fs << RED[0] << RED[1] << RED[2];
					break;
				case(3) :
					fs << GREEN[0] << GREEN[1] << GREEN[2];
					break;

				default:
					fs << WHITE[0] << WHITE[1] << WHITE[2];
					break;
				}
			}
		}
		fs.close();
	}

	void traceplot(string filename, mat val) {
		int len = val.size();
		double mean = 0;
		double max = val(0,0);
		double min = val(0,0);
		for (int i = 0; i < len; i++) {
			mean += val(i, 0);
			if (max < val(i, 0)) max = val(i, 0);
			if (min > val(i, 0)) min = val(i, 0);
		}
		mean = mean / len;
		double range = max - min;
		if (range < .00000000001) range = 1;
		vector<int> xval(len);
		vector<int> yval(len);
		for (int i = 0; i < len; i++) {
			int u = floor(double(i) / double(len) * 500);
			xval[i] = u;
			yval[i] = floor( 500*(val(i, 0)-min)/range);
		}
		BITMATWRITE500 trplot(filename, xval, yval);
		trplot.BITMAPCREATING_CT();
	}
}
