/*
Copyright (C) 2012 David Doria, daviddoria@gmail.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "ImageGraphCut.h"

// Submodules
#include "Mask/ITKHelpers/Helpers/Helpers.h"
#include "Mask/ITKHelpers/ITKHelpers.h"
#include "Mask/StrokeMask.h"

// ITK
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkTestingComparisonImageFilter.h"
#include "itkVectorImage.h"

int main(int argc, char*argv[])
{
  // Verify arguments
  if(argc != 5)
  {
    std::cerr << "Required: image.png foreground.stroke background.stroke baseline.png" << std::endl;
    return EXIT_FAILURE;
  }

  // Parse arguments
  std::string imageFilename = argv[1];

  std::string foregroundStrokeFilename = argv[2];

  std::string backgroundStrokeFilename = argv[3];

  std::string baselineFilename = argv[4];

  // Output arguments
  std::cout << "imageFilename: " << imageFilename << std::endl
            << "foregroundStrokeFilename: " << foregroundStrokeFilename << std::endl
            << "backgroundStrokeFilename: " << backgroundStrokeFilename << std::endl
            << "baselineFilename: " << baselineFilename << std::endl;

  // The type of the image to segment
  typedef itk::VectorImage<float, 2> ImageType;

  // Read the image
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(imageFilename);
  reader->Update();

  // Read the foreground and background stroke images
  StrokeMask::Pointer foregroundStrokeMask =
      StrokeMask::New();
  foregroundStrokeMask->Read(foregroundStrokeFilename);

  StrokeMask::Pointer backgroundStrokeMask =
      StrokeMask::New();
  backgroundStrokeMask->Read(backgroundStrokeFilename);

  // Perform the cut
  ImageGraphCut<ImageType> GraphCut;
  GraphCut.SetImage(reader->GetOutput());
  GraphCut.SetNumberOfHistogramBins(20);
  GraphCut.SetLambda(.01);
  std::vector<itk::Index<2> > foregroundPixels =
      ITKHelpers::GetPixelsWithValue(foregroundStrokeMask.GetPointer(), StrokeMaskPixelTypeEnum::STROKE);
  std::cout << "There are " << foregroundPixels.size() << " foreground pixels." << std::endl;

  std::vector<itk::Index<2> > backgroundPixels =
      ITKHelpers::GetPixelsWithValue(backgroundStrokeMask.GetPointer(), StrokeMaskPixelTypeEnum::STROKE);
  std::cout << "There are " << backgroundPixels.size() << " background pixels." << std::endl;

  GraphCut.SetSources(foregroundPixels);
  GraphCut.SetSinks(backgroundPixels);
  GraphCut.PerformSegmentation();

  // Get and write the result
  ForegroundBackgroundSegmentMask* result = GraphCut.GetSegmentMask();

  // Read the baseline image
  ForegroundBackgroundSegmentMask::Pointer baselineMask =
      ForegroundBackgroundSegmentMask::New();
  baselineMask->Read(baselineFilename);

  baselineMask->Write("BaselineReadWrite.png", ForegroundPixelValueWrapper<unsigned char>(255),
                BackgroundPixelValueWrapper<unsigned char>(0));

  unsigned int numberOfIncorrectPixels = ITKHelpers::CountDifferentPixels(result, baselineMask.GetPointer());
  if(numberOfIncorrectPixels > 0)
  {
    std::cerr << "There were " << numberOfIncorrectPixels << " incorrect pixels!" << std::endl;
    result->Write("Incorrect.png", ForegroundPixelValueWrapper<unsigned char>(255),
                  BackgroundPixelValueWrapper<unsigned char>(0));
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
