
import Tkinter

w = None

def on_enter_item(e):
    print 'ohai', w.find_withtag(Tkinter.CURRENT)

def on_leave_item(e):
    print 'kbye', w.find_withtag(Tkinter.CURRENT)

def on_mouse_release(e):
    print 'mouse release', w.canvasx(e.x), w.canvasy(e.y)

def on_mouse_press(e):
    print 'mouse press', w.canvasx(e.x), w.canvasy(e.y)

def on_mouse_motion(e):
    print 'mouse motion', w.canvasx(e.x), w.canvasy(e.y)

def on_delete(e):
    print 'delete'

def main():

    global w

    master = Tkinter.Tk()

    w = Tkinter.Canvas(master, width=640, height=480, highlightthickness=0)
    w.pack(fill=Tkinter.BOTH, expand=1)

    hover_color = 'snow2'
    normal_color = 'snow3'

    hover_color_conflict = 'coral2'
    normal_color_conflict = 'coral3'

    hover_color_merge = 'skyblue2'
    normal_color_merge = 'skyblue3'

    #dash_pattern = [3,3]
    dx, dy = (50, 25)
    element_a = w.create_rectangle(50, 25, 150, 75, fill=normal_color, activefill=hover_color)
    element_b = w.create_rectangle(50 + dx, 25 + dy, 150 + dx, 75 + dy, fill=normal_color_conflict, activefill=hover_color_conflict)
    element_c = w.create_rectangle(50 + 2*dx, 25 + 2*dy, 150 + 2*dx, 75 + 2*dy, fill=normal_color_merge, activefill=hover_color_merge)

    w.bind('<ButtonPress-1>', on_mouse_press)
    w.bind('<ButtonRelease-1>', on_mouse_release)
    w.bind('<B1-Motion>', on_mouse_motion)
    w.bind('<Delete>', on_delete)

    for el in (element_a, element_b, element_c):
        w.tag_bind(el, '<Any-Enter>', on_enter_item)
        w.tag_bind(el, '<Any-Leave>', on_leave_item)

    w.focus_set()

    Tkinter.mainloop()

if __name__ == '__main__':
    main()
