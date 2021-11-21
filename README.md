## Install and run project
```
git clone https://github.com/vubom01/smoke_simulation
cd smoke_simulation
mkdir build
cd build
cmake ..
make
./smoke_simulator
```
## Cấu trúc project
```
├── CGL                     // https://github.com/Computer-graphics-INT3403/hw_0/tree/main/CGL
├── ext/nanogui             // https://github.com/wjakob/nanogui
├── src
│   ├── callback.cpp        // xử lý thao tác với chuột và bàn phím
│   ├── callback.h          // header 
│   ├── color.h             // chuyển đổi màu sắc từ hsv sang rgb và ngược lại
│   ├── common.cpp          // cấu hình chung project
│   ├── common.h            // header 
│   ├── grid.cpp            // khởi tạo và xử lý mô phỏng khói trên lưới   
│   ├── grid.h              // header 
│   ├── main.cpp            // main
│   ├── shader.h            // xử lý bóng
│   ├── smoke_screen.cpp    // bảng màu và thuộc tính
│   └── smoke_screen.h      // header  
├── .gitignore  
├── build.sh
├── CMakeLists.txt
└── README.md  
```
## Keyboard
| Key  | Function         |
| ---- | -----------------|
| p    | pause/resume     |
| r    | reset all params |