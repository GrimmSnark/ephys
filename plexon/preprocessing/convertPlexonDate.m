function plexonStructure = convertPlexonDate(plexonStructure)

% plexonStructure.DateTime = sprintf('%4d%2d%2d%2d%2d%2d',  plexonStructure.CreatorDateTime.Year+1900, plexonStructure.CreatorDateTime.Month+1, plexonStructure.CreatorDateTime.MonthDay, plexonStructure.CreatorDateTime.Hour, plexonStructure.CreatorDateTime.Minute, plexonStructure.CreatorDateTime.Second);
plexonStructure.DateTime = str2num(sprintf('%4d%2d%2d%02d%2d%2d',  plexonStructure.CreatorDateTime.Year+1900, plexonStructure.CreatorDateTime.Month+1, plexonStructure.CreatorDateTime.MonthDay, plexonStructure.CreatorDateTime.Hour, plexonStructure.CreatorDateTime.Minute, plexonStructure.CreatorDateTime.Second));

end