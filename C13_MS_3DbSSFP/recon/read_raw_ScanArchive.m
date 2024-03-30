function rawdata = read_raw_ScanArchive(path)

tmp = dir([path 'Scan*.h5']); 
archive = GERecon('Archive.Load', [tmp(end).folder '/' tmp(end).name]);
ControlCountAll = archive.ControlCount;
NumTR = ControlCountAll - 1;
for n = 1:ControlCountAll
    currentControl = GERecon('Archive.Next',archive);
    if currentControl.opcode > 0
        rawdata(:,:,n) = double(currentControl.Data);
    end
end
GERecon('Archive.Close', archive);

end