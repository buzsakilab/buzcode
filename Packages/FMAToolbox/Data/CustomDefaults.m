%CustomDefaults - User-defined custom default values for function properties.
%
%  The behavior of most FMAToolbox functions can be modified using optional
%  property-value pairs, e.g. GetPositions returns position samples in different
%  coordinate systems depending on the property 'coordinates'. Such properties
%  usually have default values. In some cases, default values can be customized.
%
%  Note that customizing a default value is *different* from assigning an explicit
%  value to a property: it applies when the function is called without explicitly
%  setting a value for the property.
%
%  An example will make this clear. By default, GetPositions returns position
%  samples in normalized coordinates. To use video coordinates instead, one can
%  of course call:
%
%    p = GetPositions('coordinates','video');
%
%  but this does not affect the *default* behavior of GetPositions, which is
%  to return normalized coordinates if one calls:
%
%    p = GetPositions;
%
%  What custom default values allow is to override this default behavior and
%  have GetPositions return video coordinates when called without explicit
%  parameters.
%
%  To define a custom default value, users can e.g. add the following lines to
%  their startup.m file:
%
%   global SETTINGS;
%   SETTINGS.GetPositions.video = 'video';
%
