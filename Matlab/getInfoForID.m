function [height,weight,age,gender,pants] = getInfoForID(fileName, ID)
    D = strsplit(ID,'-');
    ID = D{1};
    C = readtable(fileName);
    for i=1:length(C.Var1)
        if(C.Var1{i} == ID)
            height = C.Var6(i);
            weight = C.Var3(i);
            age = C.Var5(i);
            gender = C.Var4(i);
            pants = C.Var9(i);
            return;
        end
    end
end