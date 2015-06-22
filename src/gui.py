#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
#
# Author:  Mario S. Könz <mskoenz@gmx.net>
# Date:    23.06.2014 10:36:54 CEST
# File:    gui.py

# I only documented the non obviouse parts of the gui, since one needs to be familiar with QT in order to work on this file. 

from qt_import import *
from .integration_routine import integrate_cf
from .plot import plot

from src_import import *

class Q2LaserWidget(QWidget):
    __laser_counter = 0
    def __init__(self, parent = None):
        self.__class__.__laser_counter += 1
        super(Q2LaserWidget, self).__init__(parent)
        self.init_ui()
    def close(self):
        self.__class__.__laser_counter -= 1
    def init_ui(self):
        #=================== widgets ===================
        self.name = QLineEdit(self)
        self.freq = QLineEdit(self)
        self.phase = QLineEdit(self)
        self.pattern = QLineEdit(self)
        #------------------- settings ------------------- 
        self.name.setPlaceholderText("name")
        self.freq.setPlaceholderText("freq")
        self.phase.setPlaceholderText("phase")
        self.pattern.setPlaceholderText("pattern")
        
        self.name.setStatusTip("name of the laser")
        self.freq.setStatusTip("frequency of the laser")
        self.phase.setStatusTip("phase of the laser")
        self.pattern.setStatusTip("firing pattern of the laser e.g: [[0, pi], [2*pi, 4.1]]")
        
        self.pattern.setMinimumWidth(130)
        
        self.name.setText("laser" + str(self.__class__.__laser_counter))
        #=================== layout ===================
        grid = QGridLayout()
        grid.addWidget(self.name   , 1, 1, 1, 1)
        grid.addWidget(self.freq   , 1, 2, 1, 1)
        grid.addWidget(self.phase  , 1, 3, 1, 1)
        grid.addWidget(self.pattern, 1, 4, 1, 1)
        
        self.setLayout(grid)
        
        self.show()

class Q2AtomWidget(QWidget):
    __atom_counter = 0
    def __init__(self, parent = None):
        self.__class__.__atom_counter += 1
        super(Q2AtomWidget, self).__init__(parent)
        self.init_ui()
    def close(self):
        self.__class__.__atom_counter -= 1
    
    def init_ui(self):
        #=================== widgets ===================
        self.name = QLineEdit(self)
        self.levels = QLineEdit(self)
        self.state = QLineEdit(self)
        
        #------------------- settings ------------------- 
        self.name.setPlaceholderText("name")
        self.levels.setPlaceholderText("energy levels")
        self.state.setPlaceholderText("start state")
        
        self.name.setStatusTip("name of the atom")
        self.levels.setStatusTip("energy levels e.g: [0, 10, 20]")
        self.state.setStatusTip("start state e.g: diag([1, 0, 0])")
        self.levels.setMinimumWidth(100)
        self.state.setMinimumWidth(100)
        
        self.name.setText("atom" + str(self.__class__.__atom_counter))
        #=================== layout ===================
        grid = QGridLayout()
        grid.addWidget(self.name   , 1, 1, 1, 1)
        grid.addWidget(self.levels , 2, 1, 1, 1)
        grid.addWidget(self.state  , 3, 1, 1, 1)
        
        self.setLayout(grid)
        self.show()
    
class Q2VibronWidget(QWidget):
    __vibron_counter = 0
    def __init__(self, parent = None):
        self.__class__.__vibron_counter += 1
        super(Q2VibronWidget, self).__init__(parent)
        self.init_ui()
    def close(self):
        self.__class__.__vibron_counter -= 1
    
    def init_ui(self):
        #=================== widgets ===================
        self.name = QLineEdit(self)
        self.freq = QLineEdit(self)
        self.state = QLineEdit(self)
        
        #------------------- settings ------------------- 
        self.name.setPlaceholderText("name")
        self.freq.setPlaceholderText("freq")
        self.state.setPlaceholderText("start state")
        
        self.name.setStatusTip("name of the vibron")
        self.freq.setStatusTip("frequency of the vibron")
        self.state.setStatusTip("start state e.g: diag([1, 0, 0])")
        self.state.setMinimumWidth(130)
        
        self.name.setText("vib" + str(self.__class__.__vibron_counter))
        #=================== layout ===================
        grid = QGridLayout()
        grid.addWidget(self.name   , 1, 1, 1, 1)
        grid.addWidget(self.freq , 1, 2, 1, 1)
        grid.addWidget(self.state  , 1, 3, 1, 1)
        
        self.setLayout(grid)
        
        self.show()

class Q2IntegrationWidget(QWidget):
    def __init__(self, parent = None):
        super(Q2IntegrationWidget, self).__init__(parent)
        
        #=================== widgets ===================
        self.hbar = QLineEdit(self)
        self.interval = QLineEdit(self)
        self.measure = QLineEdit(self)
        self.cutoff = QLineEdit(self)
        self.adv_phi = QLineEdit(self)
        
        #------------------- settings ------------------- 
        self.hbar.setPlaceholderText("hbar")
        self.hbar.setText("1.05457173e-34")
        self.interval.setPlaceholderText("intervall")
        self.measure.setPlaceholderText("N")
        self.cutoff.setPlaceholderText("cutoff")
        self.adv_phi.setPlaceholderText("advanced phi")
        
        self.hbar.setStatusTip("value for the reduced planck constant")
        self.interval.setStatusTip("integration intervall e.g: [0, 10*pi]")
        self.measure.setStatusTip("numbers of measurements")
        self.cutoff.setStatusTip("cutoff frequency for rwa")
        self.adv_phi.setStatusTip("input an entangled phi")
        
        #=================== layout ===================
        grid = QGridLayout()
        
        grid.addWidget(self.hbar         , 0, 0, 1, 1)
        grid.addWidget(self.interval     , 0, 1, 1, 1)
        grid.addWidget(self.measure      , 0, 2, 1, 1)
        grid.addWidget(self.cutoff       , 0, 3, 1, 1)
        grid.addWidget(self.adv_phi      , 0, 4, 1, 1)
        self.setLayout(grid)
        
        self.show()

class Q2NumberWidget(QWidget):
    def __init__(self, parent = None):
        super(Q2NumberWidget, self).__init__(parent)
        
        #=================== widgets ===================
        self.atom_sb = QSpinBox(self)
        self.laser_sb = QSpinBox(self)
        self.vibron_sb = QSpinBox(self)
        
        self.atom_sb.setRange(1, 10)
        self.laser_sb.setRange(1, 20)
        self.vibron_sb.setRange(0, 10)
        #------------------- settings ------------------- 
        self.atom_sb.setStatusTip("number of atoms")
        self.laser_sb.setStatusTip("number of lasers")
        self.vibron_sb.setStatusTip("number of vibrons")
        
        #=================== layout ===================
        grid = QGridLayout()
        
        grid.addWidget(QLabel("Laser") , 1, 1, 1, 1)
        grid.addWidget(QLabel("Atom")  , 2, 1, 1, 1)
        grid.addWidget(QLabel("Vibron"), 3, 1, 1, 1)
        grid.addWidget(self.laser_sb  , 1, 2, 1, 1)
        grid.addWidget(self.atom_sb   , 2, 2, 1, 1)
        grid.addWidget(self.vibron_sb , 3, 2, 1, 1)
        self.setLayout(grid)
        
        self.show()

class Q2PlotWidget(QWidget):
    def __init__(self, parent = None):
        super(Q2PlotWidget, self).__init__(parent)
        grid = QGridLayout(self)
        
        self.parent = parent
        self.text = QTextEdit(self)
        self.btn = QPushButton("Delete and Autogenerate", self)
        self.btn.pressed.connect(self.autogen)
        
        grid.addWidget(self.text, 1, 1, 1, 1)
        grid.addWidget(self.btn, 2, 1, 1, 1)
        self.setLayout(grid)
        self.show()
    
    def autogen(self):
        text = self.parent.autogen()
        self.text.setText(text)
    
class Q2DisplayWidget(QMainWindow):
    def __init__(self, parent = None):
        super(Q2DisplayWidget, self).__init__(parent)
        self.laser_w = []
        self.atom_w = []
        self.vibron_w = []
        self.rabi_w = []
        self.eta_w = []
        self.init_ui()
        
    def init_ui(self):
        #=================== widgets ===================
        self.itg_w = Q2IntegrationWidget(self)
        self.num_w = Q2NumberWidget(self)
        self.plot = Q2PlotWidget(self)
        self.collapse = QTextEdit(self)
        
        #=================== connects ===================
        self.num_w.laser_sb.valueChanged.connect(self.changed_sb)
        self.num_w.atom_sb.valueChanged.connect(self.changed_sb)
        self.num_w.vibron_sb.valueChanged.connect(self.changed_sb)
        
        #=================== actions ===================
        exit_action = QAction(self.style().standardIcon(QStyle.SP_TitleBarCloseButton), "Exit", self)
        exit_action.setShortcut("Ctrl+Q")
        exit_action.triggered.connect(self.close)

        load_action = QAction(self.style().standardIcon(QStyle.SP_DialogOpenButton), "Open", self)
        load_action.setShortcut("Ctrl+O")
        load_action.triggered.connect(self.load)
        
        save_action = QAction(self.style().standardIcon(QStyle.SP_DialogSaveButton), "Save as", self)
        save_action.setShortcut("Ctrl+S")
        save_action.triggered.connect(self.save)
        
        itg_action = QAction(self.style().standardIcon(QStyle.SP_MediaPlay), "Integrate", self)
        itg_action.setShortcut("Ctrl+P")
        itg_action.triggered.connect(self.integrate)
        
        #=================== menu bar ===================
        self.menuBar()
        menubar = self.menuBar()
        file_menu = menubar.addMenu("&File")
        file_menu.addAction(load_action)
        file_menu.addAction(save_action)
        file_menu.addAction(itg_action)
        file_menu.addAction(exit_action)
        
        #=================== tool bar ===================
        self.toolbar = self.addToolBar('File')
        #~ self.toolbar.addAction(exit_action)
        self.toolbar.addAction(load_action)
        self.toolbar.addAction(save_action)
        self.toolbar.addAction(itg_action)
        
        #=================== layout ===================
        self.init_dynaminc_ui()
        
        self.statusBar()
        self.setWindowTitle("working title")
    
    def init_dynaminc_ui(self):
        self.setup = QWidget(self)
        grid = QGridLayout()
        grid.addWidget(self.num_w     , 2, 1, 1, 1)
        
        #------------------- create lasers ------------------- 
        while self.num_w.laser_sb.value() < len(self.laser_w):
            x = self.laser_w.pop()
            x.close()
            del x
            y = self.rabi_w.pop()
            del y
        while self.num_w.laser_sb.value() > len(self.laser_w):
            self.laser_w.append(Q2LaserWidget(self))
            self.rabi_w.append([])
            for i in self.atom_w:
                self.rabi_w[-1].append(QLineEdit())
                self.rabi_w[-1][-1].setStatusTip("rabi frequency (pass a list e.g. [1, 2] if the rabi-freq are different for different transitions)")
        
        for l, l_i in zip(self.laser_w, range(len(self.laser_w))):
            grid.addWidget(l, 3 + l_i, 1, 1, 1)
        
        #------------------- create atoms ------------------- 
        while self.num_w.atom_sb.value() < len(self.atom_w):
            x = self.atom_w.pop()
            x.close()
            del x
            for r_row in self.rabi_w:
                y = r_row.pop()
                del y
            for e_row in self.eta_w:
                y = e_row.pop()
                del y
        while self.num_w.atom_sb.value() > len(self.atom_w):
            self.atom_w.append(Q2AtomWidget(self))
            for r_row in self.rabi_w:
                r_row.append(QLineEdit())
                r_row[-1].setStatusTip("rabi frequency (pass a list e.g. [1, 2] if the rabi-freq are different for different transitions)")
            for e_row in self.eta_w:
                e_row.append(QLineEdit())
                e_row[-1].setStatusTip("lamb-dicke parameter")
        
        for a, a_i in zip(self.atom_w, range(len(self.atom_w))):
            grid.addWidget(a, 2, 2 + a_i, 1, 1)
        
        #------------------- create vibrons ------------------- 
        while self.num_w.vibron_sb.value() < len(self.vibron_w):
            x = self.vibron_w.pop()
            x.close()
            del x
            y = self.eta_w.pop()
            del y
        while self.num_w.vibron_sb.value() > len(self.vibron_w):
            self.vibron_w.append(Q2VibronWidget(self))
            self.eta_w.append([])
            for i in self.atom_w:
                self.eta_w[-1].append(QLineEdit())
                self.eta_w[-1][-1].setStatusTip("lamb-dicke parameter")
        for v, v_i in zip(self.vibron_w, range(len(self.vibron_w))):
            grid.addWidget(v, 3 + v_i + len(self.laser_w), 1, 1, 1)
        
        #------------------- create rabi boxes ------------------- 
        for r_row, r_row_i in zip(self.rabi_w, range(len(self.rabi_w))):
            for r, r_i in zip(r_row, range(len(r_row))):
                grid.addWidget(r, 3 + r_row_i, 2 + r_i, 1, 1)
        #------------------- create eta boxes ------------------- 
        for e_row, e_row_i in zip(self.eta_w, range(len(self.eta_w))):
            for e, e_i in zip(e_row, range(len(e_row))):
                grid.addWidget(e, 3 + len(self.laser_w) + e_row_i, 2 + e_i, 1, 1)
        
        grid.addWidget(self.itg_w     , 3 + len(self.laser_w) + len(self.vibron_w), 1, 1, 2)
        self.setup.setLayout(grid)
        
        tab = QTabWidget(self)
        tab.addTab(self.setup, "Setup")
        tab.addTab(self.plot, "Plot")
        tab.addTab(self.collapse, "Collapse")
        self.setCentralWidget(tab)
        self.resize(0, 0)
        self.show()
        
    def changed_sb(self, val):
        self.init_dynaminc_ui()
    
    def autogen(self):
        res = []
        cf = self.parse()
        for a in cf["atom"]:
            res.append([])
            for lvl, i_lvl in zipi(a[1]):
                name = a[0]
                res[-1].append([name + "(fock_dm(N, " + str(i_lvl) + "))", name + " " + str(lvl)])
        
        for a in cf["vibron"]:
            res.append([])
            for i_lvl in range(len(a[3])):
                name = a[0]
                res[-1].append([name + "(fock_dm(N, " + str(i_lvl) + "))", name + " " + str(i_lvl)])
        
        res = "\n\n".join(["\n".join(["; ".join(x) for x in y]) for y in res])
        return res
    
    def parse(self):
        # parses the gui content into the canonical form for later use of the integrate and plot fct
        from .useable_fct import *
        
        cf = {}
        #=================== parse lasers ===================
        cf["laser"] = []
        for l in self.laser_w:
            name = str(l.name.text())
            freq = eval(str(l.freq.text()))
            phase = eval(str(l.phase.text()))
            pattern = eval(str(l.pattern.text()))
            cf["laser"].append([name
                              , freq
                              , phase
                              , pattern
                               ])
        #=================== parse atoms ===================
        cf["atom"] = []
        for a in self.atom_w:
            name = str(a.name.text())
            levels = eval(str(a.levels.text()))
            state = eval(str(a.state.text()))
            if isinstance(state, q.Qobj):
                state = state.data.todense()
            assert len(state) == len(levels), 'dimension mismatch in ' + name
            cf["atom"].append([name
                             , levels
                             , state
                              ])
        #=================== parse vibrons ===================
        cf["vibron"] = []
        for v in self.vibron_w:
            name = str(v.name.text())
            freq = eval(str(v.freq.text()))
            state = eval(str(v.state.text()))
            if isinstance(state, q.Qobj):
                state = state.data.todense()
            cf["vibron"].append([name
                               , freq
                               , len(state)
                               , state
                               ])
        
        #=================== parse rabi ===================
        cf["rabi"] = []
        for r_row, r_row_i in zip(self.rabi_w, range(len(self.rabi_w))):
            temp = []
            for r, r_i in zip(r_row, range(len(r_row))):
                rabi = eval(str(r.text()))
                n_lvl = len(cf["atom"][r_i][1])
                trans = n_lvl * (n_lvl - 1) / 2
                if isinstance(rabi, list):
                    assert len(rabi) == trans, 'dimension mismatch in rabi freq ' + str(rabi)
                else:
                    rabi = [rabi] * trans
                temp.append(rabi)
            
            cf["rabi"].append(temp)
        
        #=================== parse eta ===================
        cf["eta"] = []
        for i in self.atom_w:
            cf["eta"].append([])
        for e_row, e_row_i in zip(self.eta_w, range(len(self.eta_w))):
            for e, e_i in zip(e_row, range(len(e_row))):
                eta = eval(str(e.text()))
                cf["eta"][e_i].append(eta)
        
        #=================== parse integration ===================
        cf["hbar"] = eval(str(self.itg_w.hbar.text()))
        cf["upper"] = eval(str(self.itg_w.interval.text()))[1]
        cf["lower"] = eval(str(self.itg_w.interval.text()))[0]
        cf["measure"] = eval(str(self.itg_w.measure.text()))
        cf["cutoff"] = eval(str(self.itg_w.cutoff.text()))
        phi = str(self.itg_w.adv_phi.text())
        if phi == "":
            cf["adv_phi"] = ""
        else:
            cf["adv_phi"] = eval(phi)
        
        #=================== parse plot ===================
        pl = str(self.plot.text.toPlainText())
        
        pl = pl.split("\n")
        
        new = [[]]
        for p in pl:
            if p.find(";") == -1:
                if len(new[-1]) > 0:
                    new.append([])
            else:
                new[-1].append(p.split("; "))
        
        if new[-1] == []:
            new.pop()
        
        cf["plot"] = new
        
        #=================== parse collapse ===================
        col = str(self.collapse.toPlainText())
        col = col.split("\n")
        
        new = []
        for c in col:
            if c.find("(") != -1:
                new.append(c.split("#")[0])
        
        cf["collapse"] = new
        
        return cf
    
    def save_parse(self):
        # only saves the text for saving, since we want pi to stay pi and not change to 3.14 if loaded again
        # the notation str(widget.text()) instead of just widget.text() is necessary bc of PyQt4. .text() returns a str in PySide but a QString in PyQt4, therefore the cast
        cf = {}
        #=================== parse lasers ===================
        cf["laser"] = []
        for l in self.laser_w:
            name = str(l.name.text())
            freq = str(l.freq.text())
            phase = str(l.phase.text())
            pattern = str(l.pattern.text())
            cf["laser"].append([name
                              , freq
                              , phase
                              , pattern
                               ])
        #=================== parse atoms ===================
        cf["atom"] = []
        for a in self.atom_w:
            name = str(a.name.text())
            levels = str(a.levels.text())
            state = str(a.state.text())
            cf["atom"].append([name
                             , levels
                             , state
                              ])
        #=================== parse vibrons ===================
        cf["vibron"] = []
        for v in self.vibron_w:
            name = str(v.name.text())
            freq = str(v.freq.text())
            state = str(v.state.text())
            cf["vibron"].append([name
                               , freq
                               , state
                               ])
        
        #=================== parse rabi ===================
        cf["rabi"] = []
        for r_row in self.rabi_w:
            temp = []
            for r in r_row:
                rabi = str(r.text())
                temp.append(rabi)
            
            cf["rabi"].append(temp)
        
        #=================== parse eta ===================
        cf["eta"] = []
        for i in self.atom_w:
            cf["eta"].append([])
        for e_row in self.eta_w:
            for e, e_i in zip(e_row, range(len(e_row))):
                eta = str(e.text())
                cf["eta"][e_i].append(eta)
        
        #=================== parse integration ===================
        cf["interval"] = str(self.itg_w.interval.text())
        cf["measure"] = str(self.itg_w.measure.text())
        cf["cutoff"] = str(self.itg_w.cutoff.text())
        cf["adv_phi"] = str(self.itg_w.adv_phi.text())
        
        #=================== parse plot/hbar/collapse ===================
        cf["plot"] = str(self.plot.text.toPlainText())
        cf["hbar"] = str(self.itg_w.hbar.text())
        cf["collapse"] = str(self.collapse.toPlainText())
        return cf
        
    def load_parse(self, cf):
        # just set the text to the right widgets
        self.num_w.laser_sb.setValue(len(cf["laser"]))
        self.num_w.atom_sb.setValue(len(cf["atom"]))
        self.num_w.vibron_sb.setValue(len(cf["vibron"]))
        #=================== parse lasers ===================
        for l_w, l in zip(self.laser_w, cf["laser"]):
            l_w.name.setText(l[0])
            l_w.freq.setText(l[1])
            l_w.phase.setText(l[2])
            l_w.pattern.setText(l[3])
        #=================== parse atoms ===================
        for a_w, a in zip(self.atom_w, cf["atom"]):
            a_w.name.setText(a[0])
            a_w.levels.setText(a[1])
            a_w.state.setText(a[2])
        #=================== parse vibrons ===================
        for v_w, v in zip(self.vibron_w, cf["vibron"]):
            v_w.name.setText(v[0])
            v_w.freq.setText(v[1])
            v_w.state.setText(v[2])
            
        #=================== parse rabi ===================
        for r_w_row, r_row in zip(self.rabi_w, cf["rabi"]):
            temp = []
            for r_w, r in zip(r_w_row, r_row):
                r_w.setText(r)
        
        #=================== parse eta ===================
        for e_w_row, e_w_row_i in zip(self.eta_w, range(len(self.eta_w))):
            for e_w, e_w_i in zip(e_w_row, range(len(e_w_row))):
                e_w.setText(cf["eta"][e_w_i][e_w_row_i])
        
        #=================== parse integration ===================
        self.itg_w.interval.setText(cf["interval"])
        self.itg_w.measure.setText(cf["measure"])
        self.itg_w.cutoff.setText(cf["cutoff"])
        self.itg_w.adv_phi.setText(cf["adv_phi"])
        
        #=================== parse plot/hbar/collapse ===================
        self.plot.text.setText(cf["plot"])
        self.itg_w.hbar.setText(cf["hbar"])
        self.collapse.setText(cf["collapse"])
    
    def save(self):
        # note the small difference between PySide and PyQt4
        print("save")
        fileName = QFileDialog.getSaveFileName(self, "Save File", "~/", "Pickle Files (*.pickle)")
        
        if qt_binding == "PySide":
            file_name = unicode(fileName[0])
        else:
            file_name = unicode(fileName)
        
        if file_name == "":
            return
            
        if file_name.find(".pickle") == -1:
            file_name += ".pickle"
        
        data = self.save_parse()
        f = open(file_name, 'wb')
        pickle.dump(data, f)
        f.close()
        
        self.setWindowTitle(file_name.split("/")[-1])
    
    def load(self):
        # note the small difference between PySide and PyQt4
        print("load")
        fileName = QFileDialog.getOpenFileName(self, "Open File", "~/", "Pickle Files (*.pickle)")
        if qt_binding == "PySide":
            file_name = unicode(fileName[0])
        else:
            file_name = unicode(fileName)
            
        if file_name == "":
            return
            
        f = open(file_name, 'rb')
        data = pickle.load(f)
        f.close()
        data["adv_phi"] = data.get("adv_phi", "")
        self.load_parse(data)
        
        self.setWindowTitle(file_name.split("/")[-1])
    
    def integrate(self):
        
        canonical_form = self.parse()
        
        times, exp = integrate_cf(canonical_form)
        plot(times, exp, canonical_form)
    
def run_gui():
    app = QApplication(sys.argv)
    w = Q2DisplayWidget()
    app.exec_()
