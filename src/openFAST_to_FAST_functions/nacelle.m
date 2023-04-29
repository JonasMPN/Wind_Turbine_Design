function FAST_object = nacelle(FAST_object, shaft_tilt, hub_radius)
    FAST_object.Nacelle.Hub.ShaftTilt = shaft_tilt;
    FAST_object.Nacelle.Housing.Diameter = 2*hub_radius;
end

