#ifndef IMAGE_H
#define IMAGE_H

#include "basic.hpp"

#include <algorithm>
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <cmath>
#include <cfloat>


// Class to represent images, capable of reading/writing in PPM and BMP formats.
// All external coordinates will be expressed with the origin at the bottom
// left, even though most image encodings start at the top left.
class Image 
{
protected:
    int width, height;
	double *color;

public:
    // Constructor which does nothing; for calling readPPM afterwards.
    Image() : width(0), height(0), color(NULL) {}

    // Constructor which initializes the data buffer.
    Image(int width, int height) : width(width), height(height) 
	{ 
		color = new double[3 * width * height];
	}

    // Destructor.
	~Image() { delete[] color; }

	void writeBMP(std::string const &filename) const
	{
		unsigned char bmpfileheader[14] = { 'B', 'M', 0, 0, 0, 0, 0, 0, 0, 0, 54, 0, 0, 0 };
		unsigned char bmpinfoheader[40] = { 40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 24, 0 };
		//unsigned char bmppad[3] = { 0, 0, 0 };
		int filesize = 54 + 3 * width*height;

		bmpfileheader[2] = (unsigned char)(filesize);
		bmpfileheader[3] = (unsigned char)(filesize >> 8);
		bmpfileheader[4] = (unsigned char)(filesize >> 16);
		bmpfileheader[5] = (unsigned char)(filesize >> 24);

		bmpinfoheader[4] = (unsigned char)(width);
		bmpinfoheader[5] = (unsigned char)(width >> 8);
		bmpinfoheader[6] = (unsigned char)(width >> 16);
		bmpinfoheader[7] = (unsigned char)(width >> 24);
		bmpinfoheader[8] = (unsigned char)(height);
		bmpinfoheader[9] = (unsigned char)(height >> 8);
		bmpinfoheader[10] = (unsigned char)(height >> 16);
		bmpinfoheader[11] = (unsigned char)(height >> 24);

		double maxIntensity = DBL_MIN, minIntensity = DBL_MAX;
		for (int i = 0; i < 3 * width*height; i++)
		{
			if (maxIntensity < color[i])
			{
				maxIntensity = color[i];
			}
			if (minIntensity > color[i])
			{
				minIntensity = color[i];
			}
		}

		FILE *f = NULL;
		f = fopen(filename.c_str(), "wb");
		if (!f) { puts("can't write output image to disk!"); return; }

		fwrite(bmpfileheader, 1, 14, f);
		fwrite(bmpinfoheader, 1, 40, f);
		unsigned char offsetBuf = 0;
		long lineOffset = width * 3 % 4;
		if (lineOffset != 0)
		{
			lineOffset = 4 - lineOffset;
		}
		unsigned char mappedIntensity;
		for (int i = height - 1; i >= 0; i--)	//bmp stores images upside down
		{
			for (int j = 0; j<width; j++)
			{
				int start = 3 * (i*width + j);
				
				//!!! maybe a better tone mapping algorithm
				mappedIntensity = (unsigned char)((color[start + 2] - minIntensity) / (maxIntensity - minIntensity) * 255);
				//mappedIntensity = (unsigned char)(max(0.0, min(1.0, color[start + 2] / 16)) * 255);
				fwrite(&mappedIntensity, 1, 1, f);

				mappedIntensity = (unsigned char)((color[start + 1] - minIntensity) / (maxIntensity - minIntensity) * 255);
				//mappedIntensity = (unsigned char)(max(0.0, min(1.0, color[start + 1] / 16)) * 255);
				fwrite(&mappedIntensity, 1, 1, f);

				mappedIntensity = (unsigned char)((color[start] - minIntensity) / (maxIntensity - minIntensity) * 255);
				//mappedIntensity = (unsigned char)(max(0.0, min(1.0, color[start] / 16)) * 255);
				fwrite(&mappedIntensity, 1, 1, f);
			}

			fwrite(&offsetBuf, 1, lineOffset, f);
		}

		fclose(f);
	}

    // Tell us if the Image is ready to be used.
    bool good() const { return color != NULL; }

    // Look up the color (as an RGB vector) at the given UV coordinates (from 0 to 1).
	Vector lookup(double u, double v) const
	{
		// Clamp coords from 0 to 1.
		u = std::min(1.0, std::max(0.0, u));
		v = std::min(1.0, std::max(0.0, v));

		// Calculate pixel coords.
		int x = static_cast<int>(u * (width - 1));
		int y = static_cast<int>(height - v * (height - 1) - 1);

		// Calculate offset into data buffer.
		double *p = color + 3 * (y * width + x);

		// Convert to doubles from 0 to 1.
		return Vector(p[0], p[1], p[2]);
	}

    // Set a given pixel. Values should be from 0 to 1.
	void setPixel(int x, int y, Vector const& color)
	{
		y = height - y - 1;
		int offset = 3 * (y * width + x);
		for (int i = 0; i < 3; i++) {
			this->color[offset + i] = color[i];
		}
	}
	void setPixel(int x, int y, double gray) { setPixel(x, y, Vector(gray, gray, gray)); }

	int getWidth(void) { return width; }
	int getHeight(void) { return height; }
	double *getData() { return color; }
};

#endif
