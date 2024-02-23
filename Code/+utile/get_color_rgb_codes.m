function [colors] = get_color_rgb_codes(imagename)
colors = cell(1, length(imagename));

%Match input label to color output

for i =1:length(imagename)

    if strcmp('Battlefield', imagename{i})
        colors{i} = [255, 176, 0]/255;
    elseif strcmp('Swimming',  imagename{i}) 
        colors{i} = utile.rgb('Silver');    

    elseif strcmp('Spoon',  imagename{i})
        colors{i} = [254, 97, 0]/255;

    elseif strcmp('Nifzig',  imagename{i}) 
       colors{i} = utile.rgb('green');

    elseif strcmp('Cowboy',  imagename{i})  
        colors{i} = [220, 38, 127]/255;

    elseif strcmp('Python',  imagename{i}) 
        colors{i} = utile.rgb('Silver');    
        colors{i} = [120, 94, 240]/255;

    elseif strcmp('Telephone',  imagename{i}) 
        colors{i} = utile.rgb('Purple');
    elseif strcmp('Bindip',  imagename{i}) 
        colors{i} = utile.rgb('LightGreen');

    elseif strcmp('ImaginedSpeech',  imagename{i})
        colors{i} = [34,122,165]/255; 

    elseif strcmp('Speech',  imagename{i}) ==1 
        colors{i} = [30,136,229]/255; %[18,61,132]/255;
        
    elseif strcmp('ImaginedOverlapSpeech',  imagename{i})  
        colors{i} = [157,154,153]/255;

    elseif strcmp('Auditory',  imagename{i}) 
        colors{i} = [30,136,229]/255;

    elseif  strcmp('Written',  imagename{i})
        colors{i} = [93,191,59]/255;

    elseif strcmp('ITI',  imagename{i}) 
        colors{i} = [154,155,156]/255;

    elseif strcmp('Cue',  imagename{i}) 
         colors{i} = [92,179,255]/255;

    elseif strcmp('Delay1',  imagename{i}) 
        colors{i} = [220, 38, 127]/255;

    elseif strcmp('Internal',  imagename{i}) 
        colors{i} = [34,122,165]/255;

    elseif strcmp('Delay2',  imagename{i})
        colors{i} = [220, 38, 127]/255;

   elseif strcmp('Speech',  imagename{i}) 
        colors{i} =  [18,61,132]/255;
    elseif strcmp('DecodedWord',  imagename{i}) 
        colors{i} =  [18,61,132]/255;
    elseif strcmp('SMG',  imagename{i}) 
        colors{i} = [122,144,183]/255;

    elseif strcmp('S1X',  imagename{i}) 
        colors{i} = [211,65,92]/255;
    else  
        error(['Condition ' imagename{i} ' not present, add color'])
    end 
end 

end

