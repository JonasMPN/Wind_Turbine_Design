function FAST_object = tower(FAST_object, tower_data, hub_height, eff_density)
    FAST_object.Tower.Height = tower_data.height;
    FAST_object.Tower.Diameter = tower_data.diameter;
    FAST_object.Tower.WallThickness = tower_data.wall_thickness;
    FAST_object.Tower.HubHeight = hub_height;
    FAST_object.Tower.EffectiveDensity = eff_density;
end

