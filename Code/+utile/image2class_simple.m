function [ class_all ] = image2class_simple( imagename_all )

%Returns the name of a class if provided a number, or the number of the
%class if provided a name. 

if ischar(imagename_all)
    image_n = 1;
else
    image_n = length(imagename_all); 
end

class_all = {};

for n_words = 1:image_n

    if ~isnumeric(imagename_all) 

        if image_n ~= 1
            imagename = imagename_all{n_words};
        else
            imagename = imagename_all;
        end 

        if strcmp('Battlefield', imagename)
             class =1;

         elseif strcmp('Cowboy', imagename)
             class =2;

        elseif strcmp('Python', imagename)
             class =3;

        elseif strcmp('Spoon', imagename)
              class =4;

        elseif strcmp('Swimming', imagename)
               class =5;

        elseif strcmp('Telephone', imagename)
               class =6;

        elseif strcmp('Bindip', imagename)
               class =7;

        elseif strcmp('Nifzig', imagename)
                class =8;
        else

            error([ imagename ' - Unknown label, add it to list']);
        end
        class_all{n_words} = class; 
        
        
    elseif(isnumeric(imagename_all))
        imagename = imagename_all(n_words);

        if imagename ==1 
             class = 'Battlefield';

        elseif imagename == 2
             class ='Cowboy';

        elseif  imagename == 3
             class = 'Python';

        elseif imagename ==4
              class ='Spoon'; 

        elseif imagename== 5
               class = 'Swimming';

        elseif imagename ==6
               class ='Telephone';

        elseif imagename == 7
               class ='Bindip';

        elseif imagename == 8
                class ='Nifzig';                   
        else
            error([ imagename 'Unknown grasp, add it to list']);
        end
        class_all{n_words} = class; 

    end 

end 


if isnumeric(class_all{1})
    class_all = cell2mat(class_all);
end




end

