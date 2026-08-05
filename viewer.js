class OrbitalViewer extends HTMLElement {
  constructor() {
    super();
    this.electrons = [];
    this.selectedElectron = this.selectedSlice = 0;
    this.innerHTML = `<div class="button-row"></div><div class="energy-label"></div><div class="slider-row"><input id="slice" type="range"><span id="slice-label"></span></div><canvas id="grid" width="420" height="420"></canvas><div id="hover">Hover over a cell to inspect its value</div><p id="hint"></p>`;
  }

  set data(electrons) {
    this.electrons = electrons;
    this.selectedElectron = this.selectedSlice = 0;
    this.render();
  }

  connectedCallback() {
    this.canvas = this.querySelector("#grid");
    this.querySelector("#slice").oninput = (ev) => {
      this.selectedSlice = +ev.target.value;
      this.updateSliceLabel();
      this.draw();
    };
    this.canvas.onmousemove = (ev) => this.onHover(ev);
    this.canvas.onmouseleave = () => this.querySelector("#hover").textContent = "Hover over a cell to inspect its value";
  }

  render() {
    if (!this.electrons.length) return;
    const e = this.electrons[this.selectedElectron];
    const [nx, ny, nz] = e.dims;
    const btn = this.querySelector(".button-row");
    btn.innerHTML = this.electrons.map((el, i) =>
      `<button class="electron-btn ${i === this.selectedElectron ? "selected" : ""}" data-i="${i}">Electron ${i + 1}${el.spin ? ` (${el.spin})` : ""}</button>`
    ).join("");
    btn.querySelectorAll("button").forEach((b) => b.onclick = () => {
      this.selectedElectron = +b.dataset.i;
      this.render();
    });
    this.querySelector(".energy-label").textContent = `Energy: ${e.energyHartree.toFixed(6)} Hartree`;
    const sl = this.querySelector("#slice");
    sl.max = nz - 1;
    sl.value = this.selectedSlice;
    this.updateSliceLabel();
    this.querySelector("#hint").textContent = `Grid: ${nx} x ${ny} x ${nz} points, spacing ${e.spacingNm.toFixed(4)} nm`;
    this.draw();
  }

  updateSliceLabel() {
    const e = this.electrons[this.selectedElectron];
    const offset = (this.selectedSlice - (e.dims[2] - 1) / 2) * e.spacingNm;
    this.querySelector("#slice-label").textContent = `Slice ${this.selectedSlice + 1} / ${e.dims[2]} (z ≈ ${offset.toFixed(3)} nm)`;
  }

  draw() {
    const e = this.electrons[this.selectedElectron];
    const [nx, ny] = e.dims;
    const max = Math.max(...e.density.flat(2)) || 1e-12;
    const ctx = this.canvas.getContext("2d");
    const cw = 420 / nx, ch = 420 / ny;
    for (let x = 0; x < nx; x++) {
      for (let y = 0; y < ny; y++) {
        const t = Math.min(1, e.density[x][y][this.selectedSlice] / max);
        ctx.fillStyle = `rgb(${Math.round(255 * Math.min(1, t * 2))},${Math.round(255 * (1 - Math.abs(t - 0.5) * 2))},${Math.round(255 * Math.min(1, (1 - t) * 2))})`;
        ctx.fillRect(x * cw, y * ch, cw, ch);
      }
    }
  }

  onHover(ev) {
    const e = this.electrons[this.selectedElectron];
    const [nx, ny] = e.dims;
    const rect = this.canvas.getBoundingClientRect();
    const x = Math.floor((ev.clientX - rect.left) / 420 * nx);
    const y = Math.floor((ev.clientY - rect.top) / 420 * ny);
    if (x >= 0 && x < nx && y >= 0 && y < ny) {
      const v = e.density[x][y][this.selectedSlice];
      this.querySelector("#hover").textContent = `x=${x}, y=${y}, z=${this.selectedSlice} → ${v.toExponential(6)}`;
    }
  }
}

customElements.define("orbital-viewer", OrbitalViewer);