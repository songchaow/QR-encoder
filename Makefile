qrcode: main.c definition.h
	clang main.c -I /home/songchaow/Codes/allegro_linux/allegro-5.0.10/installation/include -L/home/songchaow/Codes/allegro_linux/allegro-5.0.10/installation/lib -lm -lallegro-static -lallegro_primitives-static -lpthread -lX11 -lGL -lXcursor -o qrcode
qrcodestt: main.c definition.h
	clang main.c -I /home/songchaow/Codes/allegro_linux/allegro-5.0.10/installation/include -L/home/songchaow/Codes/allegro_linux/allegro-5.0.10/installation/lib -lm -lallegro-static -lallegro_primitives-static -lpthread -l:libX11.a -lGL -l:libXcursor.a -v -o qrcode