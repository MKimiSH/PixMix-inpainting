function [res] = inpaint_video_PixMix(invid, maskvid)
% invid: t x h x w x ch
% mask: t x h x w
res = zeros(size(invid));

prevFr = [];
currFr = [];

prevMap = [];
currMap = [];

for fr = 1:size(invid, 1)
   currFr = squeeze(invid(fr, :,:,:));
   currMask = squeeze(maskvid(fr, :,:));
   if isempty(prevFr) || isempty(prevMap)
      [currFr, currMap] = fillOneImage(currFr, mask, 0);
   elseif tooDifferent(currFr, prevFr)
      [currFr, currMap] = fillOneImage(currFr, mask, 0);
   else
      [currFr, currMap] = fillOneImage_withinit(currFr, mask, prevFr, prevMap)
   end
   
   prevFr = currFr;
end

end