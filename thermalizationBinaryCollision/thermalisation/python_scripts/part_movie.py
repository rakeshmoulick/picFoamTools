import cv2
import os

def save():
    os.system("ffmpeg -r 2 -i ./figure_part/part_%d.png -vcodec mpeg4 -y movie_part.mp4") # the no. shows the speed
    #os.system("ffmpeg -r 10 -i ./figure_den/part_%d.png -vcodec mpeg4 -y movie_den.mp4")

save()