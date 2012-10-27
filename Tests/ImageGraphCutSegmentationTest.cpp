#include "ImageGraphCut.h"

// Submodules
#include "Mask/ITKHelpers/Helpers/Helpers.h"
#include "Mask/ITKHelpers/ITKHelpers.h"

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkVectorImage.h"

int main(int argc, char*argv[])
{
  // Verify arguments
  if(argc != 5)
    {
    std::cerr << "Required: image foregroundMask backgroundMask output" << std::endl;
    return EXIT_FAILURE;
    }

  // Parse arguments
  std::string imageFilename = argv[1];

  // This image should have white pixels indicating foreground pixels and be black elsewhere.
  std::string foregroundFilename = argv[2];

  // This image should have white pixels indicating background pixels and be black elsewhere.
  std::string backgroundFilename = argv[3];

  std::string outputFilename = argv[4]; // Foreground/background segment mask

  // Output arguments
  std::cout << "imageFilename: " << imageFilename << std::endl
            << "foregroundFilename: " << foregroundFilename << std::endl
            << "backgroundFilename: " << backgroundFilename << std::endl
            << "outputFilename: " << outputFilename << std::endl;

  // The type of the image to segment
  typedef itk::VectorImage<float, 2> ImageType;

  // Read the image
  typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(imageFilename);
  reader->Update();

  // Read the foreground and background stroke images
  Mask::Pointer foregroundMask = Mask::New();
  foregroundMask->Read(foregroundFilename);

  Mask::Pointer backgroundMask = Mask::New();
  backgroundMask->Read(backgroundFilename);

  // Perform the cut
  ImageGraphCut<ImageType> GraphCut;
  GraphCut.SetImage(reader->GetOutput());
  GraphCut.SetNumberOfHistogramBins(20);
  GraphCut.SetLambda(.01);
  std::vector<itk::Index<2> > foregroundPixels = ITKHelpers::GetNonZeroPixels(foregroundMask.GetPointer());
  std::vector<itk::Index<2> > backgroundPixels = ITKHelpers::GetNonZeroPixels(backgroundMask.GetPointer());
  GraphCut.SetSources(foregroundPixels);
  GraphCut.SetSinks(backgroundPixels);
  GraphCut.PerformSegmentation();

  // Get and write the result
  Mask* result = GraphCut.GetSegmentMask();

  ITKHelpers::WriteImage(result, outputFilename);
}
