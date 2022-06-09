# ITK JCW tutorial - part 1
Notebook and code (binder links) for the tutorial is on:
https://github.com/vicory/jcmsk

Contents:
- [Introduction to ITK](https://github.com/KitwareMedical/2019-03-13-KRSCourseInBiomedicalImageAnalysisAndVisualization/blob/6d55ff7ebf0f79ce62d1f8d0ba9547f9273c0d50//1_Introduction_to_the_Insight_Toolkit.ipynb) [15’ talk]
- [ITK in Python](https://github.com/KitwareMedical/2019-03-13-KRSCourseInBiomedicalImageAnalysisAndVisualization/blob/6d55ff7ebf0f79ce62d1f8d0ba9547f9273c0d50//4_ITK_in_Python.ipynb) [75’ syntax-oriented notebook]

### Function calls
- The common way to call and pipeline an ITK function is:
```python
# smooth an image
smoothed = itk.curvature_flow_image_filter(image,
                                          number_of_iterations=6,
                                          time_step=0.005)
view(smoothed, slicing_planes=True, gradient_opacity=0.8, ui_collapsed=True)
```

