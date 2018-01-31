function [circleContainer] = Plot3D_ElecBound_01( eleStruct )

% OLD plot3DboundaryElec

zINDS = cellfun(@(x) ~isempty(x), eleStruct.meanCirSMint(:,1));
circlesAll = eleStruct.meanCirSMint(zINDS);
circlesAll = cellfun(@(x) x(:,1:3), circlesAll, 'UniformOutput', false);

circleContainer = [];

for i = 1:length(circlesAll)
   
    tmp = circlesAll{i};
    
    circleContainer = [circleContainer ; tmp];

end

boundaryObj = boundary(circleContainer, 0.5);
hold on

trisurf(boundaryObj,circleContainer(:,1),circleContainer(:,2),circleContainer(:,3),'Facecolor','black','FaceAlpha',0.3,'Edgecolor','none')

xlim([0 512])
ylim([0 512])

end

