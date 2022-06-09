# ITK JCW tutorial - part 1
Notebook and code (binder links) for the tutorial is on:
https://github.com/vicory/jcmsk

Contents:
- [Introduction to ITK](https://github.com/KitwareMedical/2019-03-13-KRSCourseInBiomedicalImageAnalysisAndVisualization/blob/6d55ff7ebf0f79ce62d1f8d0ba9547f9273c0d50//1_Introduction_to_the_Insight_Toolkit.ipynb) [15’ talk]
- [ITK in Python](https://github.com/KitwareMedical/2019-03-13-KRSCourseInBiomedicalImageAnalysisAndVisualization/blob/6d55ff7ebf0f79ce62d1f8d0ba9547f9273c0d50//4_ITK_in_Python.ipynb) [75’ syntax-oriented notebook]

#### Function calls
- The common way to call and pipeline an ITK function is:
```python
# smooth an image
smoothed = itk.curvature_flow_image_filter(image,
                                          number_of_iterations=6,
                                          time_step=0.005)
view(smoothed, slicing_planes=True, gradient_opacity=0.8, ui_collapsed=True)
```

#### To print the supported types, run the following command in your python environment:
```python
itk.CurvatureFlowImageFilter.GetTypes()
```

#### Printing image metadata
```python
print(myimage)
```

#### The functions to convert ITK images to/from NumPy arrays/views are:
```python
itk.array_view_from_image()
```
and
```python
itk.image_view_from_array()
```
You can see the keyword view in both the names of these functions. ITK will not create a physical copy of the image.
Removing the keywork view will create local copies of the image.
```python
image_view = itk.GetImageViewFromArray( np_array)
image = itk.GetImageFromArray( np_array)

array_view = itk.GetArrayViewFromImage(image)
array = itk.GetArrayFromImage(image)
```
**Warning: 3D orientation is different between ITK (, simpleITK) and ndarrays! You will need to transpose the data.**
#### ITK Image and NumPy array index order
Let's look at the size of an ITK image and the size of a NumPy array.
```python
im=itk.imread("data/CBCT-TextureInput.png", itk.UC)
arr = itk.array_view_from_image(im)
print("Image size: " + str(itk.size(im)))
print("Array size: " + str(arr.shape))
```
```markdown
Image size: itkSize2 ([570, 326])
Array size: (326, 570)
```
The sizes appear to be inverted!

#### Cast to a different image type:
```python
InputImageType = itk.Image[itk.F, 2]
OutputImageType = itk.Image[itk.UC, 2]
cast_filter = itk.CastImageFilter[InputImageType, OutputImageType].New(image)
cast_filter.Update()

itk.CurvatureFlowImageFilter(cast_filter.GetOutput())
```

### Create your own pipeline
```python
def my_func(ImageType):
    my_pipeline = itk.pipeline()
    mean_filter = itk.MeanImageFilter[ImageType, ImageType].New()
    my_pipeline.connect(mean_filter)
    my_pipeline.expose("Radius")
    threshold_filter = itk.ThresholdImageFilter[ImageType].New()
    my_pipeline.connect(threshold_filter)
    my_pipeline.expose("Lower")
    return my_pipeline

PixelType = itk.UC
ImageType = itk.Image[PixelType, 2]
image = itk.imread("data/CBCT-TextureInput.png", PixelType)

my_pipeline = my_func(ImageType)
my_pipeline.SetInput(image)
my_pipeline.SetRadius(5)
my_pipeline.SetLower(100)
my_pipeline.Update()

view(my_pipeline.GetOutput())
```