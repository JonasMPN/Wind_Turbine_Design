function FAST_object = nacelle(FAST_object, shaft_tilt, hub_diameter)
    FAST_object.Nacelle.Hub.ShaftTilt = shaft_tilt;
    FAST_object.Nacelle.Housing.Diameter = hub_diameter;
end

